"""
redocking_validation.py
=========================
Validacao do protocolo de docking por REDOCKING do ligante co-cristalizado.

O QUE ISSO SIGNIFICA E POR QUE IMPORTA
-----------------------------------------
Cada estrutura do PDB usada como alvo (ex.: 1EVE, a AChE) foi resolvida por
cristalografia JA COM um inibidor conhecido dentro do sitio ativo (o
"ligante de referencia" que usamos para definir a caixa de busca). Isso da
uma oportunidade de teste: se a gente pegar esse MESMO ligante, "esquecer"
sua posicao original, e mandar o Vina docka-lo de novo do zero na mesma
caixa -- o resultado deveria ficar bem perto de onde ele estava de verdade
no cristal. Se ficar, o protocolo de docking (preparo de receptor, caixa de
busca, parametros do Vina) e confiavel para os outros ~380 compostos que
nao tem gabarito experimental. Essa e a mesma logica de validacao usada
em estudos publicados da area (ex.: Negru et al., in vivo 2025).

O quao "perto" e medido pelo RMSD (root-mean-square deviation) entre a
pose redocada e a pose experimental original, calculado atomo-a-atomo
(pareado pela ordem dos atomos pesados, ver parse_heavy_atoms_ordered).
Convencao usada na
literatura: RMSD <= 2,0 A = alta precisao; entre 2,0 e 3,0 A = precisao
razoavel; >= 3,0 A = baixa precisao (protocolo pouco confiavel para
aquele alvo especifico).

COMO USAR
---------
py redocking_validation.py --alvos alvos\\targets_summary.csv --alvos-dir alvos --outdir validacao --vina-path vina.exe

(--alvos-dir e a mesma pasta que voce passou como --outdir para o
target_prep.py -- e onde estao os arquivos <PDBID>_raw.pdb)

SAIDAS
------
- validacao/rmsd_redocking.csv       -> RMSD e classificacao por alvo
- validacao/<PDBID>_redock.pdbqt     -> pose redocada (para inspecao/figura)
- validacao/<PDBID>_crystal_ref.pdb  -> ligante cristalografico extraido (referencia)
"""

import argparse
import subprocess
import sys
import re
from pathlib import Path

import numpy as np
import pandas as pd

AFFINITY_RE = re.compile(r"^\s*1\s+(-?\d+\.?\d*)", re.MULTILINE)


def extract_reference_ligand(raw_pdb_path: Path, ligand_resname: str, preferred_chain: str):
    """
    Mesma logica de selecao de cadeia usada no target_prep.py: se houver
    mais de uma copia do ligante (proteina multimerica), usa so uma copia
    (a da cadeia pedida, ou a primeira disponivel), para nao misturar
    2 sitios de ligacao diferentes.
    Retorna a lista de linhas HETATM daquela copia.
    """
    by_chain = {}
    for line in raw_pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line[0:6].strip() != "HETATM":
            continue
        resname = line[17:20].strip()
        if resname != ligand_resname:
            continue
        chain_id = line[21:22].strip()
        by_chain.setdefault(chain_id, []).append(line)

    if not by_chain:
        return None
    use_chain = preferred_chain if preferred_chain in by_chain else sorted(by_chain.keys())[0]
    return by_chain[use_chain]


def parse_heavy_atoms_ordered(path, only_first_model=False):
    """
    Le um PDBQT e retorna a lista (na ORDEM em que aparecem no arquivo) das
    coordenadas dos atomos PESADOS (exclui tipos de hidrogenio HD/H do
    AutoDock). Usamos ORDEM em vez de NOME porque o Open Babel renomeia os
    atomos ao converter para PDBQT (perde nomes originais como "C1", vira
    soh "C"), entao pareamento por nome nao e confiavel. Ja a ORDEM dos
    atomos pesados e preservada pelo Vina entre o ligante de entrada e a
    pose de saida (o Vina so move coordenadas, nao adiciona/remove/reordena
    atomos) -- por isso e um criterio robusto para RMSD.
    """
    coords = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if only_first_model and line.startswith("ENDMDL"):
            break
        record = line[0:6].strip()
        if record not in ("ATOM", "HETATM"):
            continue
        adtype = line[77:].strip().split()[0] if len(line) > 77 and line[77:].strip() else ""
        if adtype in ("HD", "H"):
            continue  # pula hidrogenios
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except ValueError:
            continue
        coords.append((x, y, z))
    return coords


def compute_rmsd(ref_coords, pose_coords):
    """RMSD pareado por ORDEM (heavy atoms), entre o ligante de referencia
    (coordenadas cristalograficas, antes do docking) e a pose redocada."""
    n = min(len(ref_coords), len(pose_coords))
    if n == 0:
        return None, 0
    if len(ref_coords) != len(pose_coords):
        print(f"  aviso: numero de atomos pesados diferente entre referencia ({len(ref_coords)}) "
              f"e pose redocada ({len(pose_coords)}) -- comparando so os primeiros {n}.")
    ref_xyz = np.array(ref_coords[:n])
    pose_xyz = np.array(pose_coords[:n])
    rmsd = np.sqrt(np.mean(np.sum((ref_xyz - pose_xyz) ** 2, axis=1)))
    return float(rmsd), n


def classify_rmsd(rmsd):
    if rmsd is None:
        return "indeterminado"
    if rmsd <= 2.0:
        return "alta precisao"
    if rmsd <= 3.0:
        return "precisao razoavel"
    return "baixa precisao"


def run_vina_dock(vina_path, receptor, ligand, center, size, out_pdbqt, exhaustiveness=8, timeout=300):
    cmd = [
        vina_path, "--receptor", str(receptor), "--ligand", str(ligand),
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness), "--out", str(out_pdbqt),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError:
        print("\nERRO: comando do Vina nao encontrado. Confira --vina-path.")
        sys.exit(1)
    except subprocess.TimeoutExpired:
        return None
    if result.returncode != 0:
        print(f"  erro do Vina: {result.stderr.strip()[:200]}")
        return None
    match = AFFINITY_RE.search(result.stdout)
    return float(match.group(1)) if match else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alvos", required=True, help="targets_summary.csv (do target_prep.py)")
    ap.add_argument("--alvos-dir", required=True, help="pasta onde estao os <PDBID>_raw.pdb (mesma usada no target_prep.py)")
    ap.add_argument("--outdir", default="validacao")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--vina-path", default="vina")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)
    alvos_dir = Path(args.alvos_dir)

    alvos = pd.read_csv(args.alvos)
    alvos = alvos[alvos["status"] == "ok"]

    rows = []
    for _, row in alvos.iterrows():
        pdb_id = row["pdb_id"]
        print(f"\n=== {row['alvo']} ({pdb_id}) ===")

        raw_path = alvos_dir / f"{pdb_id}_raw.pdb"
        if not raw_path.exists():
            print(f"  aviso: {raw_path} nao encontrado, pulando.")
            rows.append({"alvo": row["alvo"], "pdb_id": pdb_id, "status": "raw_pdb_nao_encontrado"})
            continue

        ref_lines = extract_reference_ligand(raw_path, row["ligand_resname"], row["chain"])
        if not ref_lines:
            print(f"  aviso: ligante de referencia nao encontrado em {raw_path}, pulando.")
            rows.append({"alvo": row["alvo"], "pdb_id": pdb_id, "status": "ligante_nao_encontrado"})
            continue

        crystal_ref_path = outdir / f"{pdb_id}_crystal_ref.pdb"
        crystal_ref_path.write_text("\n".join(ref_lines) + "\nEND\n", encoding="utf-8")

        # converte o ligante cristalografico para PDBQT (ligante flexivel, mesmo processo do ligand_prep.py)
        ligand_pdbqt = outdir / f"{pdb_id}_crystal_ref.pdbqt"
        subprocess.run(
            ["obabel", str(crystal_ref_path), "-O", str(ligand_pdbqt), "-p", "7.4", "--partialcharge", "gasteiger"],
            capture_output=True, text=True, timeout=60,
        )
        if not ligand_pdbqt.exists():
            print("  aviso: falha ao converter ligante de referencia para PDBQT.")
            rows.append({"alvo": row["alvo"], "pdb_id": pdb_id, "status": "falha_conversao_pdbqt"})
            continue

        print("  redockando o ligante de referencia na sua propria caixa de busca...")
        redock_out = outdir / f"{pdb_id}_redock.pdbqt"
        center = (row["center_x"], row["center_y"], row["center_z"])
        size = (row["size_x"], row["size_y"], row["size_z"])
        affinity = run_vina_dock(args.vina_path, row["receptor_pdbqt"], ligand_pdbqt, center, size,
                                  redock_out, args.exhaustiveness)
        if affinity is None or not redock_out.exists():
            rows.append({"alvo": row["alvo"], "pdb_id": pdb_id, "status": "falha_docking"})
            continue

        ref_atoms = parse_heavy_atoms_ordered(str(ligand_pdbqt))
        pose_atoms = parse_heavy_atoms_ordered(str(redock_out), only_first_model=True)
        rmsd, n_atomos = compute_rmsd(ref_atoms, pose_atoms)
        classificacao = classify_rmsd(rmsd)

        if rmsd is None:
            print("  aviso: nao foi possivel calcular RMSD (nenhum atomo pesado pareado).")
        else:
            print(f"  RMSD = {rmsd:.2f} A ({n_atomos} atomos pareados) -> {classificacao}")
        print(f"  afinidade da pose redocada: {affinity} kcal/mol")

        rows.append({
            "alvo": row["alvo"], "pdb_id": pdb_id, "status": "ok",
            "rmsd_angstrom": round(rmsd, 2) if rmsd is not None else None,
            "n_atomos_pareados": n_atomos, "classificacao": classificacao,
            "afinidade_redock_kcal_mol": affinity,
        })

    summary = pd.DataFrame(rows)
    summary_path = outdir / "rmsd_redocking.csv"
    summary.to_csv(summary_path, index=False)

    print(f"\n\nConcluido -> {summary_path}")
    print("\nInterpretacao: RMSD <= 2.0 A = alta precisao (protocolo validado); 2.0-3.0 A = razoavel;")
    print(">= 3.0 A = baixa precisao (revisar preparo do receptor/caixa de busca para aquele alvo).")


if __name__ == "__main__":
    main()
