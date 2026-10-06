# Data preparation pipeline:

These two dataset are new ready to be fixed about redundancy, which is removed by clustering sequences by similarity. In this way we can remove redundant sequences by keeping only the representative of each cluster. We use MMseq tool with 30% sequence-identity threshold and a 40% length-coverage threshold as hyperparameters.

After that we merged positive and negative datasets and split randomly this set in two further sets each, specifically with 80% in one called training set and the remaining (20%) in other one called test set. We used an apposite python function to make the split without generate biases with respect negative/positive labels. Then the training set is again splitted into 5 fold for cross-validation.


The entire pipeline is contained in data_collection.py, which can be run with the following command in terminal:
```
python3 data_collection.py option(pos,neg,all) output_filename
```

## Commands used for clustering:

General command: 
```
mmseqs easy-cluster input.fa cluster-results tmp --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1
```
- input.fa — input FASTA file
- cluster-results — output prefix · tmp — temporary folder
- --min-seq-id 0.3 — 30% sequence-identity threshold
- -c 0.4 — 40% length-coverage threshold
- --cov-mode 0 — coverage of query and target
- --cluster-mode 1 — connected-component clustering

Command for positive set:
```
mmseqs easy-cluster LB2_project_Group_1/Data_collection/positive/positive_set.fasta cluster-results-positive tmp_pos  --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1
```
Command for negative set: 
```
mmseqs easy-cluster ../Data_collection/negative cluster_result_negative tmp --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1
```

## DA ELIMINARE O RIVEDERE MOLTO

This set have been randomly splits into two further sets each, specifically with 80% in one called training set and the remaining (20%) in other one called test set, doing this for preserve positive/negative ratio in each subset we put ``stratify = y``, so that the training set and test set maintain the same proportion.
At this point one training and one benchmark set have been obtained, these latter ones will be used in the end for look how trained models works. In the training set is done the 5-fold CV using StratifiedKfold, that means split the set into 5 subsets, this is done for evaluating the model's performance more reliably. Its been done using a for cycle, in which with the command 
``for fold_idx, (train_index, val_index) in enumerate(skf.split(X_train, y_train), start=1)``: 
- fold_idx (for the actual number of the fold 1-5)
- train_index (list of positional indexes "ID" of the rows used for the training in this fold)
- val_index (list of positional indexes "ID" of the rows used for the validation in this fold)

``X_fold_train, X_fold_val = X_train.iloc[train_index], X_train.iloc[val_index]``
In X_fold_train were putted all the observations rows associated with the value obtained in train_index, and in X_fold_val the data associated with the value found in val_index, for the validation fold.
So at this point 80% of observations are splitted in 5 folds for the training, and the remaining 20% observations are in the fold for the validation.

Afterwards, for the dataset with only the variable "class" the same procedure has been done with:
``y_fold_train, y_fold_val = y_train.iloc[train_index], y_train.iloc[val_index]``


``actual_val_indices = X_train.index[val_index]``
In this step the positional index assigned before in the actual indexes of the original data have been converted.

In the end with the command 
``df.loc[actual_val_indices, 'fold'] = fold_idx``
the column called "fold" its been created, in which is written the number of the group to which it belongs, "-1" if it's from the test set.


