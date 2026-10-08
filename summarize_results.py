#!/usr/bin/env python3
"""Turn resultados_docking.csv into Markdown ranking tables.

Reads the CSV written by run_docking.py (columns: cid, alvo, pdb_id,
afinidade_kcal_mol, pose_pdbqt, status) and prints the top N compounds per
target, ready to paste into a README or report. Rows without a score
(failed dockings) are ignored.

Example:
    python summarize_results.py resultados_docking.csv --top 5

The column names can be changed to summarize any other results file:
    python summarize_results.py other.csv --target-col target \
        --compound-col compound --score-col affinity --top 10
"""

import argparse
import csv


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv_file")
    ap.add_argument("--target-col", default="alvo")
    ap.add_argument("--compound-col", default="cid")
    ap.add_argument("--score-col", default="afinidade_kcal_mol")
    ap.add_argument("--top", type=int, default=5, help="top N per target")
    args = ap.parse_args()

    by_target = {}
    with open(args.csv_file, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            try:
                score = float(row[args.score_col])
            except (ValueError, KeyError, TypeError):
                continue
            by_target.setdefault(row[args.target_col], []).append(
                (row[args.compound_col], score))

    for target, items in by_target.items():
        items.sort(key=lambda x: x[1])
        print(f"\n**{target}**\n")
        print("| Rank | Compound (CID) | Affinity (kcal/mol) |")
        print("|-----:|----------------|--------------------:|")
        for i, (compound, score) in enumerate(items[:args.top], start=1):
            print(f"| {i} | {compound} | {score:.2f} |")


if __name__ == "__main__":
    main()
