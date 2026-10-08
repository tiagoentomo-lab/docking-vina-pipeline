#!/usr/bin/env python3
"""Turn a docking results CSV into a Markdown table for the README.

Works with the CSV written by run_docking.py and, through the column
options, with any other results file (for example a CSV that also has
compound names or CIDs).

Example:
    python summarize_results.py results/docking_results.csv --top 5
    python summarize_results.py my_results.csv --receptor-col target \
        --ligand-col compound --score-col affinity --top 10
"""

import argparse
import csv


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv_file")
    ap.add_argument("--receptor-col", default="receptor")
    ap.add_argument("--ligand-col", default="ligand")
    ap.add_argument("--score-col", default="affinity_kcal_mol")
    ap.add_argument("--top", type=int, default=5, help="top N per receptor")
    args = ap.parse_args()

    by_receptor = {}
    with open(args.csv_file, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            try:
                score = float(row[args.score_col])
            except (ValueError, KeyError, TypeError):
                continue
            by_receptor.setdefault(row[args.receptor_col], []).append(
                (row[args.ligand_col], score))

    for receptor, items in by_receptor.items():
        items.sort(key=lambda x: x[1])
        print(f"\n**{receptor}**\n")
        print("| Rank | Ligand | Affinity (kcal/mol) |")
        print("|-----:|--------|--------------------:|")
        for i, (ligand, score) in enumerate(items[:args.top], start=1):
            print(f"| {i} | {ligand} | {score:.1f} |")


if __name__ == "__main__":
    main()
