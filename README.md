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

UniprotKB query:
```
```

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

Redundancy is removed from both dataset by using clustering. To remove nearly identical sequences, only the representative of each cluster was kept. MMseq tool was used with 30% sequence-identity threshold and a 40% length-coverage threshold as hyperparameters.

Then positive and negative dataset were merged and split randomly in trainig and test set, which respectivelly contains 80% and 20% of the entire dataset. Biases with respect the class label (positive and egative) were avoided by using specific parameter in the splitting function.


## Data Analysis:

In this step we describe the datasets with respect some features (SP lengths, protein lengths, amino-acid composition of SPs vs a background, taxonomic classification and SP cleavege sites) to ensure that no biases has been introduced during the splitting (in particular between training and test set).


Da eliminare:
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
