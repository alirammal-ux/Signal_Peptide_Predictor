# Signal Peptide Prediction in Eukaryotic Proteins

This project builds and evaluates models to predict whether a eukaryotic protein contains an N-terminal **signal peptide (SP)**, and where it is cleaved. We built the dataset from UniProtKB, analysed it, and trained two classifiers: a **Von Heijne PSWM** baseline and an **SVM** with feature extraction and selection.

## Project workflow

1. **Data collection**: positive (with SP) and negative (without SP) sets from UniProtKB.
2. **Data preparation**: redundancy removal, train/test split, cross-validation folds.
3. **Data analysis**: checks that the data are consistent and balanced across splits.
4. **Von Heijne method**: PSWM baseline and its evaluation.
5. **SVM**: feature extraction, feature selection and evaluation *(to be completed)*.

## Repository structure

| Folder | Content |
|---|---|
| `Data_collection/` | Scripts and TSV files for the positive and negative datasets (`data_collection.py`) |
| `Data_preparation/` | Clustering with MMseqs2, train/test split and 5-fold CV |
| `Data_analysis/` | Plots: amino acid composition, cleavage-site logos, protein and SP lengths, taxonomy |
| `Von_Heijne/` | PSWM implementation (`matrix.ipynb`, `von_heine_metrics.py`) and results |
| `SVM/` | SVM with feature extraction and selection *(coming soon)* |

Each folder has its own README with the details.

## Data at a glance

- **Positives**: eukaryotic, reviewed, non-fragment proteins from UniProtKB with an experimentally verified SP, length ≥ 40 aa, evidence at protein level. After filtering (SP starting at the beginning of the chain and longer than 13 aa, cleaved) we kept **2961** proteins.
- **Negatives**: eukaryotic, reviewed, non-fragment proteins with no SP, experimentally located in cytosol, nucleus, mitochondrion, plastid, peroxisome or cell membrane. Proteins of the Golgi, ER, extracellular space or annotated as secreted were excluded, leaving **19594** proteins.
- **Redundancy removal**: MMseqs2 clustering (30% sequence identity, 40% coverage), keeping one representative per cluster.
- **Splits**: stratified 80% training / 20% test, with the training set divided into 5 folds for cross-validation.

## Data analysis findings:

- Hydrophobic residues (L, A, V) dominate the SP region, and the cleavage site shows the typical **AXA motif** at positions −3/−1.
- Positive proteins are shorter than negative ones (median 274 vs 425 aa).
- SP length is about 20 residues and does not correlate with protein length.
- Amino acid composition, SP length and taxonomy (mostly Metazoa) are well balanced across the 5 folds and the test set.

## Models

### Von Heijne (PSWM)
A Position-Specific Weight Matrix built on a 15-residue window (13 positions before the cleavage site, 2 after), slid over the first 90 residues of each protein. The best window score is compared with a threshold chosen to maximise the F1 score. Results are averaged over 5 runs (MCC, F1, precision, recall, accuracy, confusion matrix) and are consistent across folds.

### SVM
*To be added.*

## How to run

Dataset collection:




## Requirements

- Python 3 (pandas, numpy, scikit-learn, matplotlib)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [WebLogo](https://weblogo.threeplusone.com/) for the sequence logos

## Authors

Group 1 – Piccoli V., Rammal A., Acazi D., Zarcone P., Giavatto L. I.