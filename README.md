# LB2_project_Group_1

## Data collection pipeline:

The preliminary dataset was collected as two set: positive dataset and negative dataset.

### Positive dataset:
For the <u>preliminary positive dataset</u>, proteins with *signal peptides* were collected from UniprotKB using the following fields:
- Signal peptide only with experimental evidence
- Reviewed
- No fragments
- Eukaryotic proteins (Taxonomy: 2759)
- Sequence length >= 40 residues
- Existence of protein at protein level

The result was 2974 entries, which were filtered again using a python script to retain just entries having a signal peptide that starts at the beginning of the chain and is strictly longer than 13 amino acids, excluding not cleaved ones. After this second filtering we obtained 2961 entries.

### Negatives dataset:

For the <u>negative dataset</u>, proteins with *signal peptides* were collected from UniprotKB using the following fields:
- Reviewed
- No fragments
- Eukaryotic proteins (Taxonomy: 2759)
- Sequence length >= 40 residues
- No signal peptides at any evidence at protein level
- Experimental located in: cytosol, nucleus, mitochondrion, plastid, perozisome, cell membrane

UniprotKB query:
```
(taxonomy_id:2759) AND (existence:1) AND (reviewed:true) AND (fragment:false) AND (length:[40 TO *]) AND ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0039) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204)) NOT (ft_signal:*) NOT ((cc_scl_term:SL-0095) OR (cc_scl_term:SL-0132) OR (cc_scl_term:SL-0112) OR (cc_scl_term:SL-0243))
```
The query explicitly exlude proteins of Golgi apparatus, endoplasmic reticulum (ER), extracellular space, or were annotated as secreted proteins.

The result was 20974 entries.

## Data collection results

| \|positive dataset\| | \| negative dataset\|
| :---: | :---: |
| 2961  | 20974 |

For each dataset we produced a tsv file with features of interest and a fasta file with Uniprot IDs and sequences. 

The tsv is organize with the following columns:
1. **CLASS**: whether the data come from the positive set or from the negative, positive ones has 1.
2. **UniProtAccession**
3. **Organism name**
4. **Eukaryotic kingdom** (Metazoa, Fungi, Plants, Other)
5. **Protein length**
6. **cleavage site position**: always zero for the negative dataset
7. **Transmembrane site**: true/false if a transmembrane helix starts within the first 90 residues
(always false for the positive dataset)



##  Data Preparation Pipeline:

These two dataset are new ready to be fixed about redundancy, this is removed before splitting for clean the entire dataset. Fristly positive and negative set have been merged creating an unique file, in which in the first column the variable class says whether the protein comes from the positive (1) or from the negative (0). For remove nearly identical sequences, idenfity representative sequences and build non-redundant tranining sets. The command used is: 
```
seqs easy-cluster input.fa cluster-results tmp \ --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1
```
That says a 30% sequence-identity threshold and a 40% length-coverage threshold were used. 


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


## Data Analysis:

In this step the aim is to describe the datasets for some features like SP lengths, protein lengths, amino-acid composition of SPs vs a background, taxonomic classification and SP cleavege sites.
In particular for proteins lengths and SP lengths histogram and density plots have been done, plus a boxplot for the first one and a scatterplot for the second one.
For the amino-acid composition of SP vs a background SwissProt a combined barplot is been done.
Pie charts and a barplot have been done for the taxonomy variable, for the kingdom and for species.
Finally for the SP cleavege sites from only the sequences containing them (positive) the cleavege site region have been aligned (-13,+2) and then a sequence logo created. 


## SP detection 
For the detection of the signal-peptide the von Heijne method it has been adopted.




In the file data_collection.py there's the entire pipeline, for use it from the terminal the usage is:
```
python3 data_collection.py option(pos,neg,all) output_filename
```
