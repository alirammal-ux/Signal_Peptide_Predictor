# Data collection pipeline:

We started collecting the preliminary dataset as two set: positive data and negative data.

## Positives:

In positive data we filtered proteins with *signal peptides* directly from UniprotKB using the following fields:
- Signal peptide only with experimental evidence
- Reviewed
- No fragments
- Eukaryotic proteins (Taxonomy: 2759)
- Sequence length >= 40 residues
- Existence of protein at protein level

and we downloaded the [URL](https://rest.uniprot.org/uniprotkb/search?compressed=false&format=json&query=%28%28ft_signal_exp%3A*%29+AND+%28reviewed%3Atrue%29+AND+%28fragment%3Afalse%29+AND+%28taxonomy_id%3A2759%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28existence%3A1%29%29&size=500) for API, this first is called <u>preliminary positive dataset</u> and are present 2974 entries, from which with the following steps we filtered for others characteristics.

We have created a function called *filter_function*: it takes in input a single entry and returns True if the entry has a signal peptide that starts at the beginning of the chain and is strictly longer than 13 amino acids, excluding not cleaved ones, and return False otherwise.

After that we have created the function called *extract_fields*, this function takes in input a single entry and return some fields of it like: UniProtAccession, Organism name, Eukaryotic kingdom (Metazoa, Fungi, Plants, Other), Protein length and the position of the signal-peptide cleavage site.

In the end we obtained 2961, in this final version we organized the tsv file as follow:
First column says whether the data come from the positive set or from the negative, positive ones has 1.
Second column is the UniProtAccession.
Third column is the Organism name.
Fourth column is the Eukaryotic kingdom (Metazoa, Fungi, Plants, Other).
Fifth column is the Protein length.
Sixth column is the position of the signal-peptide cleavage site.
Seventh column for positive is always False because it's inherent to the transmembrane filter applied for the negative dataset. 


## Negatives:

In negative data we filtered proteins with *no signal peptides* directly from UniprotKB by some fields:
- Reviewed
- No fragments
- Eukaryotic proteins (Taxonomy: 2759)
- Sequence length >= 40 residues
- No signal peptides at any evidence at protein level
- Experimental located in: cytosol, nucleus, mitochondrion, plastid, perozisome, cell membrane

and we downloaded the [URL](https://rest.uniprot.org/uniprotkb/search?format=json&query=%28%28taxonomy_id%3A2759%29+AND+%28existence%3A1%29+AND+%28reviewed%3Atrue%29+AND+%28fragment%3Afalse%29+AND+%28length%3A%5B40+TO+*%5D%29+AND+%28%28cc_scl_term_exp%3ASL-0091%29+OR+%28cc_scl_term_exp%3ASL-0191%29+OR+%28cc_scl_term_exp%3ASL-0173%29+OR+%28cc_scl_term_exp%3ASL-0039%29+OR+%28cc_scl_term_exp%3ASL-0209%29+OR+%28cc_scl_term_exp%3ASL-0204%29%29+NOT+%28ft_signal%3A*%29+NOT+%28%28cc_scl_term%3ASL-0095%29+OR+%28cc_scl_term%3ASL-0132%29+OR+%28cc_scl_term%3ASL-0112%29+OR+%28cc_scl_term%3ASL-0243%29%29%29) for API, this first is the <u>negative dataset</u> and are present 20974 entries, from which we didn't have to filter any more.

In the case of negative we had to clean better our dataset because we have got in the first case 20974 observations, but in these were included the proteins of Golgi apparatus, endoplasmic reticulum (ER), extracellular space, or were annotated as secreted proteins. So we have created a better and more complex query and from wich we obtained 1380 less extra observations.

As for positive we created a function called *extract_fields* for extract only for some fields specifically:
- UniProt accession
- Organism name
- Eukaryotic kingdom (Metazoa, Fungi, Plants, Other)
- Protein length
- Whether a transmembrane helix starts within the first 90 residues (true / false)

## Final setup of the output 
In the end we obtained exactly 19594 negative cases, in this final version we organized the tsv file as follow:
First column says whether the data come from the positive set or from the negative, negative ones has 0.
Second column is the UniProtAccession.
Third column is the Organism name.
Fourth column is the Eukaryotic kingdom (Metazoa, Fungi, Plants, Other).
Fifth column is the Protein length.
Sixth column in this negative is always 0 because is a selection feature used only for positive data.
Seventh column says wheter the transmembrane helix starts within the first 90 residues (True/False).


In the end we have created the function *get_dataset*, it takes three parameters in the case of <u>negative dataset</u>: 
+ (STRING) The URL corresponding to the specific API call
+ (FUNCTION) A function to filter entries according to specific criteria and returns directly fields we are looking for
+ (STRING) The output file name where TSV will be saved
and four parameters in the case of <u>positive dataset</u>:
+ (STRING)The URL corresponding to the specific API call
+ (FUNCTION)A filter function to filter entries according to specific criteria
+ (FUNCTION)An extract function used to get some specific fields of interest from the JSON format.
+ (STRING)The output file name where TSV will be saved



#  Data Preparation Pipeline:

These two dataset are new ready to be fixed about redundancy, this is removed before splitting for clean these set. Fristly we have merged these two and for remove nearly identical sequences, idenfity representative sequences, build non-redundant tranining sets and filter similar sequences in these. The command used is: 
"seqs easy-cluster input.fa cluster-results tmp \ --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1"
That says we have used a 30% sequence-identity threshold and a 40% length-coverage threshold. 


After that we splitted randomly this set in two further sets each, specifically with 80% in one called training set and the remaining (20%) in other one called test set.
At this point we have two training and two test set, these latter ones will be used in the end. In the two training set we have done the 5-fold CV, that means split the set into 5 subsets, this is done for evaluating the model's performance more reliably.
For preserve the positive/negative ratio in each subset AABIAMO APPLICATO UN METODOOOOOOOOOOOO


In the file data_collection.py there's the entire pipeline, for use it from the terminal the usage is:
python3 data_collection.py option(pos,neg,all) output_filename
