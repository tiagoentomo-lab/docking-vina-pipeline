# Reproducible Molecular Docking Pipeline with AutoDock Vina

A small, scriptable pipeline for virtual screening by molecular docking. It runs
a library of ligands against one or more protein targets with AutoDock Vina,
stores every pose and log, and produces a single ranked CSV of binding
affinities.

The pipeline was developed for an in silico prospecting project on
isoquinoline alkaloids, carried out with high school students at SESI Amazonas
(Brazil), which won 1st place in Natural Sciences at FICSesi 2026. The same
workflow was also used in a separate study of flavonoids against human
Topoisomerase I (PDB 1T8I).

## Why command line instead of a GUI

Earlier runs used PyRx. On Windows, its GUI failed repeatedly (AutoGrid
webservice errors, IndexError crashes). Driving Vina directly from the command
line proved more stable and, more importantly, **reproducible**: one config
file, one command, a fixed random seed, and every setting recorded in code.

## Workflow

```
PubChem (CID, SDF)  ->  ligand prep (.pdbqt)  \
                                                >  Vina docking  ->  CSV  ->  ranking
PDB structure       ->  receptor prep (.pdbqt) /
```

1. **Ligands**: structures retrieved from PubChem (REST API), converted to
   `.pdbqt` with Open Babel or AutoDockTools.
2. **Receptors**: PDB structures cleaned (water, co-factors) and converted to
   `.pdbqt`.
3. **Grid box**: centered on the crystallographic ligand position, which gives
   the most reliable results when a co-crystallized ligand exists.
4. **Docking**: `run_docking.py` runs every receptor/ligand pair.
5. **Ranking**: `summarize_results.py` builds Markdown tables of the top hits.

## Repository structure

```
.
├── run_docking.py          # batch docking, resumable, writes the results CSV
├── summarize_results.py    # CSV -> Markdown ranking tables
├── config.example.json     # receptors, grid box center and size
├── ligands/                # ligand .pdbqt files (input)
├── receptors/              # receptor .pdbqt files (input)
└── results/
    └── docking_results.csv # best affinity per receptor/ligand pair
```

## Requirements

- Python 3.8 or newer (standard library only, no extra packages)
- [AutoDock Vina](https://vina.scripps.edu/) 1.2.x (tested with 1.2.7 on Windows)
- Open Babel or AutoDockTools to prepare `.pdbqt` files

## Usage

1. Put the ligand files in `ligands/` and the receptor files in `receptors/`.
2. Copy `config.example.json` to `config.json` and fill in the receptor names,
   `.pdbqt` paths and grid box (`center` and `size`, in angstroms).
3. Run (PowerShell):

```powershell
python run_docking.py --config config.json `
    --ligands-dir ligands --out-dir results `
    --vina-path "C:\Tools\vina_1.2.7_win.exe" --cpu 4
```

On Linux or macOS, use `--vina-path vina` (or the full path) and replace the
backticks with backslashes.

Useful options:

| Option | Default | Meaning |
|--------|---------|---------|
| `--exhaustiveness` | 8 | Search thoroughness (higher is slower and more thorough) |
| `--num-modes` | 9 | Poses kept per ligand |
| `--seed` | 42 | Fixed seed, so reruns give the same scores |
| `--cpu` | 4 | CPU cores used by Vina |

The script is **resumable**: pairs that already have a valid output are
skipped, so an interrupted run can simply be started again.

To print ranking tables for this README:

```bash
python summarize_results.py results/docking_results.csv --top 5
```

## Results

Binding affinity is the best Vina score (kcal/mol). More negative means a
stronger predicted interaction.

> Replace this block with the output of `summarize_results.py` for the
> isoquinoline alkaloid project.

| Rank | Ligand | Target | Affinity (kcal/mol) |
|-----:|--------|--------|--------------------:|
| 1 | | | |
| 2 | | | |
| 3 | | | |

The complete table is in [`results/docking_results.csv`](results/docking_results.csv).

## Limitations

Docking scores are predictions, not measurements, and should be read as a way
to prioritize candidates for further study.

- **No redocking validation (RMSD) yet.** Redocking the co-crystallized ligand
  to confirm that the protocol reproduces the experimental pose is planned but
  not done.
- **DNA chains removed.** AutoDockTools silently drops non-standard residues,
  which removed the DNA chains from the Topoisomerase I structure. Results for
  that target therefore describe binding to the protein without the DNA
  context.
- **Rigid receptor.** Vina treats the receptor as rigid and uses an empirical
  scoring function, so scores are approximate.
- Next steps: RMSD validation, ADMET screening and molecular dynamics for the
  best candidates.

## References

- Trott, O. and Olson, A. J. (2010). AutoDock Vina: improving the speed and
  accuracy of docking. *Journal of Computational Chemistry*, 31(2), 455-461.
- Eberhardt, J. et al. (2021). AutoDock Vina 1.2.0: new docking methods,
  expanded force field, and Python bindings. *Journal of Chemical Information
  and Modeling*, 61(8), 3891-3898.

## Author

Nicanor Tiago Bueno Antunes, biologist (Entomology, INPA) and professor at
Escola SESI Dra. Emina Barbosa Mustafa, SESI Amazonas.
