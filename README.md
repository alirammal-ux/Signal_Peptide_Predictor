# Signal Peptide Prediction in Eukaryotic Proteins

LB2 project – Group 1. We predict whether a eukaryotic protein carries an N-terminal secretory **signal peptide (SP)**. The dataset is built from UniProtKB/Swiss-Prot, made non-redundant and split into a training and a hold-out test set; on it we compare a **von Heijne position-specific weight matrix (PSWM)** with a **support vector machine (SVM)** trained on physico-chemical features and a **convolutional neural network (CNN)** trained on one-hot encoded sequences.

## Pipeline 
Project work plan:
1. Retrive relevant datasets from UniProtKB.
2. Preprocess datasets for cross-validation and benchmarking.
3. Analyze and visualize dataset statistics.
4. Extract relevant features for classification.
5. Implement von Heijne's algorithm, SVM classifier and CNN classifier.
6. Evaluates methods via cross-validation and blind test set.
7. Discuss and report results.

```mermaid
flowchart LR
    A[DATA COLLECTION<br>UniProtKB REST API<br>filtering] --> B[DATA PREPARATION<br>MMseqs2 clustering <br>train/test split + 5 CV folds]
    B --> C[DATA ANALYSIS]
    C --> D[VALIDATION of VON HEIJNE METHOD<br>pswm]
    D --> E[Build von Heijne model]
    C --> F[VALIDATION of SVM METHOD<br>Feature extraction<br>Feature selection with RF]
    F --> G[Build svm model]
    C ---> I[CNN METHOD<br>One-hot encoding<br>Hyperparameter search with Optuna]
    E --> H[RESULTS ANALYSIS<br>error analysis]
    G --> H[RESULTS ANALYSIS<br>error analysis]
    I --> H[RESULTS ANALYSIS<br>error analysis]
```

## Repository structure

| Folder | Main files | What it does |
|---|---|---|
| [`Data_collection/`](Data_collection/) | `data_collection.py` | Downloads the positive and negative sets from UniProtKB and writes them as TSV and FASTA |
| [`Data_preparation/`](Data_preparation/) | `Data_cleaning/remove_redundancy.py`, `datasplitting.py` | Removes redundancy (MMseqs2), splits training and test set, assigns the cross-validation folds |
| [`Data_analysis/`](Data_analysis/) | one notebook per analysis | Amino acid composition, cleavage-site logos, protein and SP length, taxonomy |
| [`von_Heijne_method/`](von_Heijne_method/) | `von_heijne_metrics.py` | PSWM baseline and its cross-validated evaluation for von Heijne method |
| [`SVM_method/`](SVM_method/) | `feature_extraction.py`, `feature_selection.py` | Feature extraction, feature selection and SVM in nested cross-validation to evaluate the svm method |
| [`Models/`](Models/) | `SVM_final.py`, `von_heijne_metrics_final.py`, `sigccn-pipeline-complete.ipynb` | Train and test the final models (svm, von heijne, cnn) for benchmarching. Use cross-validation to find hyperparameters. The CNN is in a notebook. |
| [`Error_Analysis/`](Error_Analysis/) | `FP_analysis/`, `FN_analysis/`, `score_distributions/` | False positive and false negative analysis of the three models on the benchmark set |

Each folder has its own README with the details.

## Dataset

| Step | Positives | Negatives |
|---|---:|---:|
| UniProtKB query and filters | 2961 | 19594 |
| After MMseqs2 clustering (30% identity, 40% coverage) | 1104 | 8624 |
| Training set (5 cross-validation folds) | 883 | 6899 |
| Hold-out test set | 221 | 1725 |

- **Positives**: reviewed, non-fragment eukaryotic proteins of at least 40 residues with evidence at protein level and an experimentally verified SP that starts at residue 1, is at least 14 residues long and is cleaved.
- **Negatives**: same constraints, no SP, and an experimentally verified localization in cytosol, nucleus, mitochondrion, plastid, peroxisome or cell membrane. Proteins also annotated in the endoplasmic reticulum, Golgi apparatus, extracellular space or as secreted are excluded.
- **Split**: stratified 80/20 split into training and test set; the training set is divided into 5 stratified folds (seed 42), each with about 177 positives and 1380 negatives.

All the later steps read `Data_preparation/Cleaned_merged_dataset.tsv`, one row per protein:

| Column | Content |
|---|---|
| `class` | 1 = signal peptide, 0 = no signal peptide |
| `ID` | UniProt accession |
| `organism_name`, `kingdom` | organism and kingdom (Metazoa, Fungi, Plants, Other) |
| `length` | protein length |
| `cleavage_pos` | last residue of the signal peptide (0 for negatives) |
| `early_trans` | helical transmembrane segment starting before residue 90 (negatives only) |
| `fold` | 1–5 = cross-validation fold, −1 = hold-out test set |

## Methods

Final models (`Models/`): each one is trained on the whole training set and evaluated once on the hold-out test set.

**von Heijne (PSWM).** The 15 residues around the cleavage site (13 before, 2 after) of the positive training sequences give a log-odds matrix against the Swiss-Prot amino acid background, with pseudocount 1. A protein is scored by sliding the window over its first 90 residues and keeping the best score, and it is predicted as SP if the score is above the threshold T = 9.26 (average of the best-F1 thresholds of the 5 cross-validation runs).

**SVM.** Each protein is described by 37 features: 17 come from sliding-window profiles of 9 AAindex scales (hydrophobicity, charge, helix and β-sheet propensity, volume, signal-sequence helical potential, flexibility, mutability) on the first 40 residues, and 20 are the amino acid composition of the first 22 residues (the median SP length in the training set). The final model keeps the top k = 35 features by Random Forest importance and uses an RBF-kernel SVM with C = 4 and γ = 0.0143.

**CNN.** Each sequence is one-hot encoded and zero-padded (or truncated) to a fixed length of 120 residues, giving a 120×20 matrix; a binary mask marks the real residues so that padding is ignored. The network has two convolutional blocks (Conv1D, ReLU, layer normalisation, dropout) acting as motif detectors for the n-, h- and c-regions: the first is a multi-kernel block with three sets of filters of different widths, the second has filters of width 9. Each block is followed by a max-pooling of width 2 (length 120 → 60 → 30). A global max-pooling (is a motif present anywhere?) and a global average-pooling (how widespread is it?) are computed in parallel over the real residues only and concatenated into a feature vector, which goes to a hidden layer of 64 ReLU neurons with dropout and to a single output neuron; a sequence is classified as SP with a fixed threshold of 0.5. The network is trained by minimising the average cross-entropy loss with the AdamW optimiser. Final configuration (see `Models/CNN/best_params.json`): 16 filters per branch with kernel widths of 5, 9 and 15 in the first layer (48 feature maps), 64 filters in the second layer, dropout rate 0.21, learning rate 5.0 × 10⁻³, weight decay 1.0 × 10⁻⁶, batch size 32 and 12 training epochs, for a total of 45,585 trainable parameters.

## Results obtained in the validation of the methods

von Heijne: in each of the 5 runs, 3 folds build the matrix, 1 fold chooses the threshold and 1 fold is used for testing. Mean ± standard error over the 5 runs (`von_Heijne_method/k_fold_stats.csv`):

| MCC | F1 | Precision | Recall | Accuracy |
|---:|---:|---:|---:|---:|
| 0.688 ± 0.039 | 0.722 ± 0.036 | 0.708 ± 0.050 | 0.745 ± 0.088 | 0.935 ± 0.008 |

SVM: A nested cross-validation on the 5 training folds chooses the number of features *k* (top *k* by Random Forest importance) and the RBF-SVM parameters *C* and γ in the inner loop, and measures the MCC of the chosen setting in the outer loop. Mean ± standard error over the nested cross-validation (see [`SVM_method/`](SVM_method/)).

| MCC |
|---:|
| 0.855 ± 0.009 | 

## Results of final models

von Heijne (`Models/vhm/von_heijne_metrics_final.py`).
SVM (`Models/SVM/SVM_final.py`).
CNN (`Models/CNN/sigccn-pipeline-complete.ipynb`).

| Model | MCC | F1 | Precision | Recall | Accuracy |
|---:|:---:|:---:|:---:|:---:|:---:|
| Von Heijne | 0.677 | 0.714 | 0.689 | 0.742 | 0.933 (normal acc) |
| SVM | 0.894 | 0.906 | 0.917 | 0.896 | 0.943 (balanced acc) |
| CNN | 0.913 | 0.923 | 0.923 | 0.923 | 0.957 (balanced acc) |

The SVM and the CNN clearly outperform the von Heijne method on the hold-out test set, and the CNN is the best-performing model on every metric.

## Data analysis highlights

- Hydrophobic residues (L, A, V) dominate the SP region, and the cleavage site shows the typical **AXA motif** at positions −3 and −1.
- Positive proteins are shorter than negative ones (median 274 vs 425 residues).
- The median SP length is 22 residues (range 14–83) and does not correlate with protein length.
- Amino acid composition, SP length and kingdoms are balanced across the 5 folds and the test set. The two classes differ in kingdom composition: Metazoa are 79% of the positives but 52% of the negatives.

Plots and details in [`Data_analysis/`](Data_analysis/).

## Installation

```bash
git clone https://github.com/Valentina-Piccoli/LB2_project_Group_1.git
cd LB2_project_Group_1
pip install -r requirements.txt
```

The project needs Python 3.6 or later (it runs on Python 3.6.11 on the course VM). `requirements.txt` lists numpy, pandas, scikit-learn (≥ 0.24), matplotlib, seaborn, biopython and requests; to open the notebooks you also need Jupyter (`pip install notebook`). Two external tools are needed only to regenerate some intermediate files: [MMseqs2](https://github.com/soedinglab/MMseqs2) for the clustering and the [WebLogo](https://weblogo.threeplusone.com/) web server for the sequence logos.

## How to run

Scripts find their files relative to their own location, so they can be launched from any folder; the commands below are run from the repository root. Every intermediate file (including the MMseqs2 clustering output) is in the repository, so you can start from any step.

```bash
# 1. Data collection (needs internet)
python3 Data_collection/data_collection.py -all Data_collection/negative/negative_set Data_collection/positive_set

# 2. Training/test split and folds -> Data_preparation/Cleaned_merged_dataset.tsv
python3 Data_preparation/datasplitting.py Data_preparation/Data_cleaning/positive_set.clean Data_preparation/Data_cleaning/negative_set.clean

# 3. von Heijne: cross-validation, then final model
python3 von_Heijne_method/von_heijne_metrics.py --cm y --spider y
python3 Models/vhm/von_heijne_metrics_final.py

# 4. SVM: features, nested cross-validation (long run), then final model
python3 SVM_method/feature_extraction.py
python3 SVM_method/feature_selection.py
python3 Models/SVM/SVM_final.py
```

**CNN.** The CNN is in the notebook `Models/CNN/sigccn-pipeline-complete.ipynb`. It has to be run on [Kaggle](https://www.kaggle.com/) and it needs a GPU (enable it in the notebook settings).

The background amino acid frequencies used by the PSWM (`Data_analysis/aa_composition_plots/background_freq.tsv`) are computed by `aa_compositions.ipynb` from the Swiss-Prot composition statistics. The notebooks locate the repository root by themselves, so they work wherever the repository is cloned.

## Authors

Group 1 – Piccoli V., Rammal A., Acazi D., Zarcone P., Giavatto L. I.