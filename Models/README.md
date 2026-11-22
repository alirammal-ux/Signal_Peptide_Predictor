# Methods

Three predictors of N-terminal signal peptides (SP), all trained on the same training set and tested once on the hold-out benchmark set.

## von Heijne (PSWM)
- Builds a position-specific weight matrix on the window around the cleavage site (13 residues before, 2 after, 15 in total) using the positive training sequences.
- The matrix is a log-odds against the Swiss-Prot amino acid background.
- A protein is scored by sliding the window over its first 90 residues and keeping the best score.
- If the best score is above a threshold, the protein is predicted as SP. The threshold (9.26) is the average of the best-F1 thresholds found in cross-validation.

## SVM
- RBF-kernel SVM (scikit-learn `SVC`).
- 37 features from the N-terminus: 17 physicochemical descriptors of the first 40 residues (AAindex scales, e.g. hydrophobicity, charge, helix propensity) and 20 amino acid frequencies.
- Features are ranked with a Random Forest (Gini importance) and the top *k* are kept.
- *k*, *C* and γ are chosen by grid search maximising the MCC (selected: k = 35, C = 4, γ = 0.0143), evaluated with nested cross-validation.

## CNN
- Input: one-hot encoded sequence, padded or truncated to 120 residues.
- Two convolutional blocks (Conv1D, ReLU, layer normalisation, dropout), each followed by max-pooling.
- Global max- and average-pooling, then a hidden layer of 64 neurons and one output neuron (SP if output ≥ 0.5).
- Trained with cross-entropy loss and the AdamW optimiser; hyperparameters were tuned with Optuna.

See the paper for the full details, formulas and results.