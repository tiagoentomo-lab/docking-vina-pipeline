"""
run_docking.py
================
Etapa 5.5 do protocolo: docking molecular e triagem virtual multi-alvo.

O QUE ESTE SCRIPT FAZ
----------------------
Roda o AutoDock Vina para TODAS as combinacoes (composto x alvo), usando:
  - os receptores e caixas de busca gerados pelo target_prep.py
    (arquivo alvos/targets_summary.csv)
  - os ligantes preparados pelo ligand_prep.py
    (arquivo saida/ligand_prep_report.csv + pasta saida/ligands_pdbqt/)

Para cada combinacao, extrai a MELHOR pose (menor/mais negativa energia de
afinidade, em kcal/mol) e vai gravando os resultados, linha por linha, num
CSV -- essa e a "matriz composto x alvo" que alimenta as etapas 5.6 (analise
de interacoes), 5.7 (SAR) e 5.8 (modelo preditivo) do protocolo.

CHECKPOINT / RETOMAR DE ONDE PAROU
------------------------------------
Como sao ~380 compostos x 5 alvos = ~1900 dockings, isso PODE DEMORAR
BASTANTE (de minutos a algumas horas, dependendo do computador e da
--exhaustiveness escolhida). Por isso, o script:
  - grava cada resultado no CSV assim que termina (nao espera acabar tudo);
  - ao ser executado de novo, PULA automaticamente as combinacoes que ja
    estao no CSV de saida -- entao se voce precisar parar (Ctrl+C) e
    continuar depois, e so rodar o mesmo comando de novo.

RECOMENDACAO DE USO
---------------------
1. Rode primeiro um teste pequeno pra ver o tempo por docking:
   py run_docking.py --alvos alvos\\targets_summary.csv --ligantes saida\\ligand_prep_report.csv --ligantes-dir saida\\ligands_pdbqt --out resultados_docking.csv --limit-compostos 5 --exhaustiveness 8

2. Veja no console quanto tempo levou por docking (ele mostra o tempo medio
   no final) e multiplique por 380 x 5 pra estimar o tempo total. Se estiver
   demorado demais, reduza --exhaustiveness (o padrao do Vina e 8; valores
   menores, tipo 4, sao mais rapidos e ainda razoaveis para triagem inicial).

3. Quando estiver satisfeito com o tempo estimado, rode sem --limit-compostos
   para o conjunto completo. Pode deixar rodando em segundo plano.

COMO USAR
---------
py run_docking.py --alvos alvos\\targets_summary.csv --ligantes saida\\ligand_prep_report.csv --ligantes-dir saida\\ligands_pdbqt --out resultados_docking.csv

SAIDA
-----
- resultados_docking.csv  -> 1 linha por (CID, alvo), com a melhor afinidade
                              (kcal/mol) e o caminho da pose docada (pdbqt)
- poses/<alvo>/<CID>_out.pdbqt -> pose 3D do melhor resultado (para a
                              analise de interacoes na etapa 5.6)
"""

import argparse
import csv
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

AFFINITY_RE = re.compile(r"^\s*1\s+(-?\d+\.?\d*)", re.MULTILINE)
FIELDNAMES = ["cid", "alvo", "pdb_id", "afinidade_kcal_mol", "pose_pdbqt", "status"]


def run_one_docking(vina_path, receptor, ligand, center, size, out_pdbqt, exhaustiveness, timeout):
    cmd = [
        vina_path,
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--out", str(out_pdbqt),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, "timeout"
    except FileNotFoundError:
        print("\nERRO: comando 'vina' nao encontrado. Instale o AutoDock Vina e garanta "
              "que esta no PATH (mesma logica usada para o Open Babel).")
        sys.exit(1)

    if result.returncode != 0:
        return None, f"erro_vina: {result.stderr.strip()[:200]}"

    match = AFFINITY_RE.search(result.stdout)
    if not match:
        return None, "sem_afinidade_no_output"
    return float(match.group(1)), "ok"


def load_done_pairs(out_csv: Path):
    """Le o CSV de resultados ja existente (se houver) para retomar de onde parou."""
    if not out_csv.exists():
        return set()
    done = set()
    try:
        df = pd.read_csv(out_csv)
        for _, row in df.iterrows():
            done.add((str(row["cid"]), str(row["alvo"])))
    except Exception:
        pass
    return done


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alvos", required=True, help="targets_summary.csv gerado pelo target_prep.py")
    ap.add_argument("--ligantes", required=True, help="ligand_prep_report.csv gerado pelo ligand_prep.py")
    ap.add_argument("--ligantes-dir", required=True, help="pasta com os arquivos <CID>.pdbqt")
    ap.add_argument("--out", default="resultados_docking.csv", help="CSV de resultados (default: resultados_docking.csv)")
    ap.add_argument("--poses-dir", default="poses", help="pasta onde salvar as poses docadas (default: poses)")
    ap.add_argument("--exhaustiveness", type=int, default=8, help="exhaustiveness do Vina (default: 8; use 4 para triagem mais rapida)")
    ap.add_argument("--timeout", type=int, default=300, help="tempo maximo (segundos) por docking antes de desistir (default: 300)")
    ap.add_argument("--limit-compostos", type=int, default=None, help="limitar a N compostos, para teste rapido")
    ap.add_argument("--vina-path", default="vina", help="caminho do executavel do Vina, se nao estiver no PATH")
    args = ap.parse_args()

    alvos_df = pd.read_csv(args.alvos)
    alvos_df = alvos_df[alvos_df["status"] == "ok"].copy()
    if alvos_df.empty:
        print("Nenhum alvo com status 'ok' encontrado em", args.alvos)
        sys.exit(1)
    print(f"{len(alvos_df)} alvo(s) prontos para docking: {list(alvos_df['alvo'])}")

    lig_df = pd.read_csv(args.ligantes)
    lig_df = lig_df[lig_df["status"] == "ok"].copy()
    if args.limit_compostos:
        lig_df = lig_df.head(args.limit_compostos)
        print(f"Modo teste: limitando a {len(lig_df)} compostos.")
    print(f"{len(lig_df)} composto(s) prontos para docking.")

    ligantes_dir = Path(args.ligantes_dir)
    out_csv = Path(args.out)
    poses_dir = Path(args.poses_dir)
    poses_dir.mkdir(exist_ok=True)

    done_pairs = load_done_pairs(out_csv)
    if done_pairs:
        print(f"Retomando: {len(done_pairs)} combinacoes ja feitas serao puladas.")

    write_header = not out_csv.exists()
    f_out = open(out_csv, "a", newline="", encoding="utf-8")
    writer = csv.DictWriter(f_out, fieldnames=FIELDNAMES)
    if write_header:
        writer.writeheader()

    total_pairs = len(lig_df) * len(alvos_df)
    done_count = 0
    times = []
    t_start = time.time()

    for _, alvo_row in alvos_df.iterrows():
        alvo_nome = alvo_row["alvo"]
        alvo_dir = poses_dir / alvo_nome.replace("/", "_").replace(" ", "_")
        alvo_dir.mkdir(exist_ok=True)
        center = (alvo_row["center_x"], alvo_row["center_y"], alvo_row["center_z"])
        size = (alvo_row["size_x"], alvo_row["size_y"], alvo_row["size_z"])
        receptor = alvo_row["receptor_pdbqt"]

        for _, lig_row in lig_df.iterrows():
            cid = str(lig_row["cid"])
            done_count += 1

            if (cid, alvo_nome) in done_pairs:
                continue

            ligand_path = ligantes_dir / f"{cid}.pdbqt"
            if not ligand_path.exists():
                writer.writerow({"cid": cid, "alvo": alvo_nome, "pdb_id": alvo_row["pdb_id"],
                                  "afinidade_kcal_mol": "", "pose_pdbqt": "", "status": "ligante_nao_encontrado"})
                f_out.flush()
                continue

            out_pdbqt = alvo_dir / f"{cid}_out.pdbqt"
            t0 = time.time()
            affinity, status = run_one_docking(
                args.vina_path, receptor, ligand_path, center, size, out_pdbqt,
                args.exhaustiveness, args.timeout,
            )
            elapsed = time.time() - t0
            times.append(elapsed)

            writer.writerow({
                "cid": cid, "alvo": alvo_nome, "pdb_id": alvo_row["pdb_id"],
                "afinidade_kcal_mol": affinity if affinity is not None else "",
                "pose_pdbqt": str(out_pdbqt) if status == "ok" else "",
                "status": status,
            })
            f_out.flush()

            if done_count % 10 == 0 or done_count == total_pairs:
                media = sum(times) / len(times) if times else 0
                restantes = total_pairs - done_count
                eta_min = (restantes * media) / 60
                print(f"[{done_count}/{total_pairs}] {alvo_nome} / CID {cid}: "
                      f"{status} ({affinity if affinity is not None else '--'} kcal/mol) | "
                      f"media {media:.1f}s/docking | ETA ~{eta_min:.0f} min")

    f_out.close()

    total_elapsed = (time.time() - t_start) / 60
    print(f"\nConcluido em {total_elapsed:.1f} min.")
    print(f"-> {out_csv}")
    print(f"-> {poses_dir}/ (poses 3D docadas, por alvo)")
    print("\nProxima etapa: 5.6 (analise das interacoes ligante-receptor dos melhores complexos)")


if __name__ == "__main__":
    main()
