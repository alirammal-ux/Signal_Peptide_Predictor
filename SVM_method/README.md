# Feature extraction and selection (SVM)

| File | Content |
|---|---|
| `feature_extraction.py` | computes the features of every protein and writes `feature.csv` |
| `aaindex_parse.py` | reads one scale from the AAindex1 database file |
| `aaindex1`, `AAindex_scales.txt` | the AAindex1 database and the list of scales used (AAindex ID, short name, description) |
| `feature_selection.py` | feature selection and SVM in nested cross-validation, writes `results_svm_rf.txt` |

## How to run

```bash
python3 Feature_Selection/feature_extraction.py      # -> Feature_Selection/feature.csv
python3 Feature_Selection/feature_selection.py       # -> Feature_Selection/results_svm_rf.txt
```

`feature_extraction.py` accepts `--fasta` (default `von_Heijne_method/merged.fasta`), `--scale_selected_file`, `--aaindex_file` and `--out` (default `Feature_Selection/feature.csv`); all defaults are found relative to the script. `feature_selection.py` has no options. The nested cross-validation trains several thousand models and takes a long time.

## Features

All profiles are computed on the first 40 residues with a sliding window: each residue gets the mean scale value of the window centred on it (the window shrinks at the ends of the sequence and ignores unknown residues).

| Scale (AAindex ID) | Window | Features | Motivation |
|---|---:|---|---|
| Hydropathy, Kyte-Doolittle (KYTJ820101) | 5 | max; position of the max; length of the region around the max with value ≥ 2; drop from the max to residue 40 | the h-region of the SP is hydrophobic |
| Positive charge, Fauchère (FAUJ880111) | 3 | max; position of the max | the n-region is positively charged |
| Negative charge, Fauchère (FAUJ880112) | 3 | max; position of the max | the c-region contains uncharged residues |
| Normalized frequency of N-terminal helix, Chou-Fasman (CHOP780204) | 9 | max; position of the max | the h-region tends to form an α-helix |
| Normalized frequency of N-terminal β-sheet, Chou-Fasman (CHOP780208) | 9 | max; position of the max | secondary structure of the N-terminus |
| Residue volume, Bigelow (BIGC670101) | 5 | mean | the c-region contains small residues |
| Signal sequence helical potential, Argos (ARGP820102) | 9 | max; position of the max | helix propensity measured on signal sequences |
| Average flexibility, Bhaskaran-Ponnuswamy (BHAR880101) | 9 | mean | |
| Relative mutability, Dayhoff (DAYM780201) | 9 | mean | |

Positions are divided by 40, so they range from 0 to 1. These are 17 features; the other 20 are the **amino acid composition of the first 22 residues**, as frequencies (22 is the median SP length in the training set, the same in every cross-validation split; non-standard residues are not counted). In total each protein has 37 features.

## Nested cross-validation

Only the 5 training folds of `Data_preparation/Cleaned_merged_dataset.tsv` are used; the hold-out test set is not touched. `feature.csv` and the dataset are sorted by UniProt ID and checked to contain the same proteins, so that each row of the features matches its class and fold.

- **Outer loop** (5 iterations): one fold is the outer test fold, the other four are the training part.
- **Inner loop** (4 iterations on the training part): every setting of the grid is trained on three folds and evaluated on the fourth; the setting with the best mean MCC is chosen (ties go to the smallest *k*).
- A setting is a triple (*k*, *C*, γ):
  - *k*: number of features kept, from 11 to total number of features in steps of 2; they are the top *k* features by Random Forest (100 trees) Gini importance, computed on the training data of that iteration only;
  - *C* ∈ {1, 4, 10, 16};
  - γ ∈ {0.5/*k*, 1/*k*, 2/*k*, 4/*k*} for the RBF kernel.
- With the chosen setting, the features are selected again on the whole training part, the SVM is trained on it and evaluated once on the outer test fold.
- Features are standardised with the mean and standard deviation of the training data only.

`results_svm_rf.txt` contains, for each outer fold, the chosen setting, the selected features and the test MCC, plus the mean and standard deviation of the MCC.

## Results


| Test fold | mcc | hyperparams (k, C, γ) |
|:---:|:---:|:---:|
| 0 | 0.8655 | (35, 4, 0.0286) |
| 1 | 0.8430 | (23, 16, 0.0217) |
| 2 | 0.8584 | (33, 4, 0.0303) |
| 3 | 0.8625 | (27, 4, 0.0185) |
| 4 | 0.8442 | (27, 4, 0.0185) |


| | |
|:---:|:---:|
|Average mcc | 0.8547 |
| Standard deviation for mcc| 0.0094