# Data preparation

In order to pass from the two raw sets in `Data_collection/` to the non-redundant, labelled and split dataset read by all the later steps, (`Cleaned_merged_dataset.tsv`), data cleaning and data preparation steps must be done. All commands are run from the repository root.

## 1. Redundancy reduction with MMseqs2

Similar sequences in the training and in the test folds would make the performance look better than it is, so each set is clustered and only one representative per cluster is kept.

```bash
mmseqs easy-cluster Data_collection/positive_set.fasta Data_preparation/Data_cleaning/positive_set.cluster/cluster-results-positive tmp --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1
mmseqs easy-cluster Data_collection/negative/negative_set.fasta Data_preparation/Data_cleaning/negative_set.cluster/cluster_result_negative tmp --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1
```

- `--min-seq-id 0.3`: 30% sequence identity threshold
- `-c 0.4`: 40% alignment coverage threshold
- `--cov-mode 0`: coverage computed on both query and target
- `--cluster-mode 1`: connected-component clustering
- `tmp`: temporary folder

Each run writes `<prefix>_rep_seq.fasta` (the representatives), `<prefix>_cluster.tsv` (cluster membership) and `<prefix>_all_seqs.fasta` in `Data_cleaning/positive_set.cluster/` and `Data_cleaning/negative_set.cluster/`.

## 2. Keeping the representatives

```bash
python3 Data_preparation/Data_cleaning/remove_redundancy.py Data_preparation/Data_cleaning/positive_set.cluster/cluster-results-positive_rep_seq.fasta Data_collection/positive_set.tsv Data_preparation/Data_cleaning/positive_set.clean
python3 Data_preparation/Data_cleaning/remove_redundancy.py Data_preparation/Data_cleaning/negative_set.cluster/cluster_result_negative_rep_seq.fasta Data_collection/negative/negative_set.tsv Data_preparation/Data_cleaning/negative_set.clean
```

`remove_redundancy.py <rep_seq.fasta> <original.tsv> <output.clean>` copies to the output only the TSV rows whose accession (column 2) is a cluster representative. It expects the 7-column TSV written by `data_collection.py`.

| | Before | After clustering |
|---|---:|---:|
| Positives | 2961 | 1104 |
| Negatives | 19594 | 8624 |

## 3. Training/test split and cross-validation folds

```bash
python3 Data_preparation/datasplitting.py Data_preparation/Data_cleaning/positive_set.clean Data_preparation/Data_cleaning/negative_set.clean
```

`datasplitting.py` merges the two `.clean` files and:

1. splits them 80/20 into training set and hold-out test set with `train_test_split`, stratified by class so that both keep the same positive/negative ratio;
2. divides the training set into 5 folds with `StratifiedKFold` (shuffled, stratified by class);
3. stores the fold of each protein in a new `fold` column: 1–5 for the training folds, −1 for the test set.

Both steps use the random seed 42, so the split is reproducible.

| Fold | Positives | Negatives | Total |
|---|---:|---:|---:|
| 0 | 177 | 1380 | 1557 |
| 1 | 177 | 1380 | 1557 |
| 2 | 176 | 1380 | 1556 |
| 3 | 176 | 1380 | 1556 |
| 4 | 177 | 1379 | 1556 |
| Test (−1) | 221 | 1725 | 1946 |

The test set is not used by any step so far: it is kept for the final comparison of the methods.

## Output files

| File | Content |
|---|---|
| `Cleaned_merged_dataset.tsv` | one row per protein with `class`, `ID`, `organism_name`, `kingdom`, `length`, `cleavage_pos`, `early_trans`, `fold`; input of all the later steps |
| `X_train.tsv`, `y_train.tsv`, `X_test.tsv`, `y_test.tsv` | the same split with the protein table (`X_*`) and the class labels (`y_*`) in separate files; the first column is the row number in the merged table |

`DataSplitting.ipynb` runs the same steps interactively.
