# Multi-target Molecular Docking Pipeline with AutoDock Vina

Scripts from an in silico bioactivity prospecting study of isoquinoline-related
alkaloids against targets linked to Alzheimer's disease, diabetes and cancer.
The study combines virtual screening by molecular docking (AutoDock Vina),
redocking validation of the protocol, and a predictive model of docking scores.

It was carried out with high school students at Escola SESI Dra. Emina Barbosa
Mustafa (SESI Amazonas, Brazil) and won 1st place in Natural Sciences at
FICSesi 2026. The code comments and console messages are in Portuguese.

## Why command line instead of a GUI

Earlier runs used PyRx. On Windows its GUI failed repeatedly (AutoGrid
webservice errors, IndexError crashes). Driving Vina directly from the command
line proved more stable and, more importantly, reproducible: every setting is
a command-line option, results are written to CSV, and the docking run can be
stopped and resumed.

## Study workflow

```
PubChem (CIDs, 3D SDF) -> ligand prep (.pdbqt)  \
                                                 > batch docking -> CSV -> SAR, predictive model
PDB structures         -> target prep (.pdbqt)  /
                                |
              redocking validation (RMSD, per target)
```

| Stage | Script | In this repository |
|-------|--------|:------------------:|
| Compound retrieval from PubChem | `pubchem_pipeline.py` | no |
| Target preparation (receptor, grid box) | `target_prep.py` | no |
| Ligand preparation (RDKit, Open Babel) | `ligand_prep.py` | no |
| Batch docking, resumable | `run_docking.py` | yes |
| Redocking validation (RMSD) | `redocking_validation.py` | yes |
| Predictive model (Random Forest) | `predict_model.py` | no |
| SAR heatmap | `sar_heatmap.py` | no |
| Ranking tables | `summarize_results.py` | yes |

The grid box of each target is centered on the co-crystallized ligand.

## Requirements

- Python 3.8 or newer, with `pandas` and `numpy` (`pip install -r requirements.txt`)
- [AutoDock Vina](https://vina.scripps.edu/) 1.2.x (the study used 1.2.7 on Windows)
- [Open Babel](https://openbabel.org/) (`obabel` on the PATH), used by the
  redocking script to convert the reference ligand

## Input files

`run_docking.py` expects the files produced by the preparation scripts:

- `targets_summary.csv`, one row per target, with the columns `atividade`,
  `alvo`, `pdb_id`, `ligand_resname`, `chain`, `status`, `receptor_pdbqt`,
  `center_x`, `center_y`, `center_z`, `size_x`, `size_y`, `size_z`,
  `n_atomos_ligante_ref`. Only rows with `status` equal to `ok` are used. See
  [`examples/targets_summary.example.csv`](examples/targets_summary.example.csv).
- `ligand_prep_report.csv`, with at least the columns `cid` and `status` (rows
  with `status` equal to `ok` are docked). See
  [`examples/ligand_prep_report.example.csv`](examples/ligand_prep_report.example.csv).
- A folder with one `<CID>.pdbqt` file per ligand.

## Usage

Redocking validation (one run per target, using the co-crystallized ligand):

```powershell
py redocking_validation.py --alvos alvos\targets_summary.csv --alvos-dir alvos --outdir validacao --vina-path vina.exe
```

Batch docking. First a small test to measure the time per docking:

```powershell
py run_docking.py --alvos alvos\targets_summary.csv --ligantes saida\ligand_prep_report.csv --ligantes-dir saida\ligands_pdbqt --out resultados_docking.csv --limit-compostos 5 --vina-path vina.exe
```

Then the full run (without `--limit-compostos`):

```powershell
py run_docking.py --alvos alvos\targets_summary.csv --ligantes saida\ligand_prep_report.csv --ligantes-dir saida\ligands_pdbqt --out resultados_docking.csv --vina-path vina.exe
```

Each result is written to the CSV as soon as it finishes. If the run is
interrupted, running the same command again skips the pairs already done.

Useful options of `run_docking.py`:

| Option | Default | Meaning |
|--------|---------|---------|
| `--exhaustiveness` | 8 | Vina search thoroughness; lower values (for example 4) are faster for an initial screening |
| `--timeout` | 300 | Maximum seconds per docking before giving up |
| `--limit-compostos` | none | Dock only the first N compounds, for a quick test |
| `--vina-path` | `vina` | Path to the Vina executable, if it is not on the PATH |
| `--poses-dir` | `poses` | Folder for the docked poses, one subfolder per target |

Outputs:

- `resultados_docking.csv`: one row per compound and target, with the columns
  `cid`, `alvo`, `pdb_id`, `afinidade_kcal_mol` (best pose), `pose_pdbqt` and `status`
- `poses/<target>/<CID>_out.pdbqt`: docked poses
- `validacao/rmsd_redocking.csv`: RMSD and classification per target

Ranking tables from the results:

```bash
python summarize_results.py resultados_docking.csv --top 5
```

### How the redocking RMSD is computed

The co-crystallized ligand is extracted from the PDB file, converted with Open
Babel, redocked by Vina in the same grid box, and compared with its crystal
pose. Because Open Babel renames atoms during conversion, heavy atoms are
paired by their order in the file, which Vina preserves. An RMSD up to 2.0 A is
classified as high precision, up to 3.0 A as reasonable, and above that as low
precision.

## Results

Redocking validation decided which of the five targets could be interpreted:

| Target | Redocking RMSD | Decision |
|--------|---------------:|----------|
| Dipeptidyl peptidase-4 (DPP-4) | 0.23 A | Kept |
| Acetylcholinesterase (AChE) | 0.82 A | Kept |
| Tubulin | 1.12 A | Kept |
| Topoisomerase I | 10.27 A | Excluded (DNA intercalation site is incompatible with rigid-receptor docking) |
| SARS-CoV-2 main protease | not usable | Excluded (incomplete reference ligand) |

On the three validated targets, 377 of 383 compounds (98.4%, from 1,149 docking
runs) produced a valid docking, with affinities between -13.91 and -6.99
kcal/mol. The library came from 386 alkaloids catalogued in PubChem, curated to
383. More negative affinity means a stronger predicted interaction.

Best compound per validated target, and the best representative with the
isoquinoline nucleus preserved:

| Target | Compound | CID | Structural class | Affinity (kcal/mol) |
|--------|----------|----:|------------------|--------------------:|
| AChE | verticinone | 167691 | Steroidal (cevanine) | -13.91 |
| Tubulin | oxymorphone azine | 9568073 | Morphinan (dimer) | -13.82 |
| DPP-4 | cevanone glycoside | 3828935 | Steroidal (cevanine) | -11.05 |
| AChE | talicarpine | 21470 | Isoquinoline nucleus | -10.93 |
| Tubulin | coptisine | 72322 | Isoquinoline nucleus | -9.47 |
| DPP-4 | curine | 253793 | Isoquinoline nucleus | -9.35 |

Reference ligands redocked with the same protocol scored -11.00 kcal/mol
(donepezil, AChE), -9.40 kcal/mol (sitagliptin, DPP-4) and -8.20 kcal/mol
(DAMA-colchicine, tubulin). A Random Forest model with six molecular
descriptors, evaluated by 5-fold cross-validation, reached R² of 0.27 (AChE),
0.50 (tubulin) and 0.69 (DPP-4).

Other findings:

- 36 compounds (9.5%) scored better than donepezil on AChE and 25 compounds
  (6.6%) scored better than sitagliptin on DPP-4. On tubulin, 87% of the
  library scored better than DAMA-colchicine, so that criterion is weakly
  selective for this target.
- 30 compounds ranked in the top 10% for at least two targets (16 for all
  three). Of these, 26 are morphinans, mostly large dimers or semisynthetic
  conjugates, and 4 are steroidal alkaloids. Part of this multi-target profile
  probably reflects molecular size rather than real selectivity.
- Affinities differed between structural classes (Kruskal-Wallis, p < 0.001
  for all targets). Molecular weight was the most important descriptor in the
  models, and ligand efficiency analysis showed a strong size effect on the
  scores. Small, planar isoquinoline alkaloids such as coptisine are efficient
  binders for their size.

## About the compound library

Of the 377 compounds with a valid docking, 317 (84%) are morphinan derivatives,
47 (12%) have the isoquinoline nucleus, 9 (2%) are steroidal and 4 (1%) belong
to other classes. The library is therefore dominated by morphinans rather than
classic isoquinoline alkaloids, and some top hits are not isoquinolines at all
(verticinone is a steroidal alkaloid). Library curation matters when reading
the rankings.

## Limitations

Docking scores are predictions, not measurements. They help prioritize
candidates and generate hypotheses, and they do not replace in vitro or in vivo
validation.

- Vina treats the receptor as rigid and uses an empirical scoring function.
- The predictive models range from modest (AChE, R² = 0.27) to good (DPP-4,
  R² = 0.69), and part of the signal reflects molecular size.
- The library is dominated by morphinans, so conclusions do not generalize to
  all isoquinoline alkaloids.
- The model should not be extrapolated to chemical classes far from the
  training set without retraining.

A separate exploratory study of flavonoids against human Topoisomerase I (PDB
1T8I) used the same Vina workflow. In that study AutoDockTools removed the DNA
chains as non-standard residues and redocking validation was not performed, so
its results should be read with that in mind.

## References

- Trott, O. and Olson, A. J. (2010). AutoDock Vina: improving the speed and
  accuracy of docking. *Journal of Computational Chemistry*, 31(2), 455-461.
- Eberhardt, J. et al. (2021). AutoDock Vina 1.2.0: new docking methods,
  expanded force field, and Python bindings. *Journal of Chemical Information
  and Modeling*, 61(8), 3891-3898.

## Author

Nicanor Tiago Bueno Antunes, biologist (Entomology, INPA) and professor at
Escola SESI Dra. Emina Barbosa Mustafa, SESI Amazonas.
