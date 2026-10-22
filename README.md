# AHPP regression: predicting the IC50 of ACE-I inhibitory peptides

MATLAB script that predicts the IC<sub>50</sub> (µM) of peptides against angiotensin-converting enzyme-I (ACE-I) with a **regression decision tree**. It is the test model behind the Anti-hypertensive Peptide Predictor (AHPP) web server.

- Web server: <http://hazralab.iitr.ac.in/ahpp/index.php>
- Article: Kalyan et al., [*J. Agric. Food Chem.* 2021, 69, 14995–15004](https://doi.org/10.1021/acs.jafc.1c04555)

## How it works

1. **Train.** The script loads 1,466 peptides with their experimental IC<sub>50</sub> values and 31 features each. Peptides with IC<sub>50</sub> above 1000 µM are left out, leaving 1,326 for the fit, and a regression tree is trained with `fitrtree`.
2. **Read the test peptides.** Features for the new peptides come from two files, matched by peptide name.
3. **Predict.** The tree predicts the IC<sub>50</sub> of each peptide. Each is also labelled with the classes used in the paper: **positive** (≤ 100 µM), **negative** (> 1000 µM), otherwise **average**.

## Features

Each peptide has 31 features: 23 sequence-based columns followed by 8 structure-based ones.

| Group | Columns | Source |
|---|---|---|
| Sequence-based | 23 | Pseudo amino acid composition (PseAAC) from PseAAC-Builder; the paper describes 21 values (λ = 1, w = 0.05), and the files do not say what the two extra columns are |
| Structure-based | 8 | Computed from the modelled peptide structure: molecular weight, R, α, β (shape in spherical coordinates), positive charge, negative charge, hydropathy (GRAVY), volume |


## Repository contents

| Path | Purpose |
|---|---|
| [`ahpp_regression.m`](ahpp_regression.m) | Trains the tree and predicts IC<sub>50</sub> for the test peptides. |
| [`data/train/data_regr2.mat`](data/train/data_regr2.mat) | Training data: `ic50_train` (1466 × 1), `train_data1` (1466 × 23), `train_data2` (1466 × 8) and `train_data_final` (1466 × 31, the two feature blocks side by side). It also holds the peptide names and sequence tables, which MATLAB can open. |
| [`data/example_test/ml_csv.csv`](data/example_test/ml_csv.csv) | Structure features for 25 example peptides. One header line, then `sequence, mol. wt, r, alpha, beta, q_positive, q_negative, hydropathy, volume`. |
| [`data/example_test/peptides.pfpseaac.csv`](data/example_test/peptides.pfpseaac.csv) | Sequence (PseAAC) features for the same 25 peptides: 23 numbers, then the peptide name.

## Usage

Requires MATLAB with the Statistics and Machine Learning Toolbox. Open MATLAB and run:

```matlab
ahpp_regression
```

The script finds its data relative to its own location, so it can be run from any folder. It prints the predictions and writes `ic50_peptide.csv` with the columns `Peptide`, `IC50_predicted_uM` and `Class`.

To predict your own peptides, replace the two files in `data/example_test/` with files in the same format (or change the file names at the top of the script). Both files must list the same peptides; the order can differ.


## License

Copyright (C) 2026 Vivek Junghare, Hazra Lab, IIT Roorkee.

This repository is released under the **GNU General Public License v3.0**; the full text is in [`LICENSE`](LICENSE). You can use, modify and redistribute the code under those terms. Modified versions that you distribute must also be released under the GPL.
