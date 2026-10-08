#!/usr/bin/env python3
"""Batch molecular docking with AutoDock Vina.

Runs every ligand (.pdbqt) in a folder against every receptor described in a
JSON config file, keeps the Vina output and log of each pair, and writes a
single CSV with the best binding affinity (kcal/mol) per receptor/ligand pair.

The run is resumable: pairs that already have an output file are skipped, so
you can stop the script and start it again without losing work.

Example (PowerShell):
    python run_docking.py --config config.json `
        --ligands-dir ligands --out-dir results `
        --vina-path "C:\\Tools\\vina_1.2.7_win.exe" --cpu 4
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def load_config(path):
    with open(path, encoding="utf-8") as fh:
        cfg = json.load(fh)
    receptors = cfg.get("receptors", [])
    if not receptors:
        sys.exit("Config has no receptors. See config.example.json.")
    for rec in receptors:
        for key in ("name", "pdbqt", "center", "size"):
            if key not in rec:
                sys.exit(f"Receptor entry is missing '{key}': {rec}")
        if len(rec["center"]) != 3 or len(rec["size"]) != 3:
            sys.exit(f"'center' and 'size' need 3 values each: {rec['name']}")
    return receptors


def best_affinity(out_pdbqt):
    """Return the best (first) Vina score from an output .pdbqt, or None."""
    try:
        with open(out_pdbqt, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("REMARK VINA RESULT:"):
                    return float(line.split()[3])
    except (OSError, ValueError, IndexError):
        pass
    return None


def run_vina(vina, receptor, ligand, out_pdbqt, log_path, center, size, args):
    cmd = [
        str(vina),
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--exhaustiveness", str(args.exhaustiveness),
        "--num_modes", str(args.num_modes),
        "--cpu", str(args.cpu),
        "--seed", str(args.seed),
        "--out", str(out_pdbqt),
    ]
    with open(log_path, "w", encoding="utf-8") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    return proc.returncode


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config.json", help="receptor/grid box JSON")
    ap.add_argument("--ligands-dir", default="ligands", help="folder with ligand .pdbqt files")
    ap.add_argument("--out-dir", default="results", help="folder for outputs, logs and CSV")
    ap.add_argument("--vina-path", default="vina", help="path to the Vina executable")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--num-modes", type=int, default=9)
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--seed", type=int, default=42, help="fixed seed for reproducibility")
    ap.add_argument("--top", type=int, default=5, help="top N printed per receptor")
    args = ap.parse_args()

    receptors = load_config(args.config)
    ligands = sorted(Path(args.ligands_dir).glob("*.pdbqt"))
    if not ligands:
        sys.exit(f"No .pdbqt ligands found in {args.ligands_dir}")

    out_dir = Path(args.out_dir)
    poses_dir = out_dir / "poses"
    logs_dir = out_dir / "logs"
    poses_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    total = len(receptors) * len(ligands)
    done = 0
    for rec in receptors:
        for lig in ligands:
            done += 1
            tag = f"{rec['name']}__{lig.stem}"
            out_pdbqt = poses_dir / f"{tag}.pdbqt"
            log_path = logs_dir / f"{tag}.log"

            if out_pdbqt.exists() and best_affinity(out_pdbqt) is not None:
                status = "cached"
            else:
                print(f"[{done}/{total}] docking {lig.stem} -> {rec['name']}", flush=True)
                code = run_vina(args.vina_path, rec["pdbqt"], lig, out_pdbqt,
                                log_path, rec["center"], rec["size"], args)
                status = "ok" if code == 0 else f"vina_error_{code}"

            rows.append({
                "receptor": rec["name"],
                "ligand": lig.stem,
                "affinity_kcal_mol": best_affinity(out_pdbqt),
                "status": status,
            })

    csv_path = out_dir / "docking_results.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["receptor", "ligand", "affinity_kcal_mol", "status"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nResults written to {csv_path}")

    for rec in receptors:
        scored = [r for r in rows
                  if r["receptor"] == rec["name"] and r["affinity_kcal_mol"] is not None]
        scored.sort(key=lambda r: r["affinity_kcal_mol"])
        print(f"\nTop {args.top} for {rec['name']}:")
        for r in scored[:args.top]:
            print(f"  {r['ligand']:<30} {r['affinity_kcal_mol']:>7.1f} kcal/mol")
        failed = [r for r in rows if r["receptor"] == rec["name"]
                  and r["affinity_kcal_mol"] is None]
        if failed:
            print(f"  ({len(failed)} pair(s) without a score, check logs/)")


if __name__ == "__main__":
    main()
