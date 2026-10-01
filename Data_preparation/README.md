# Data preparation pipeline:

These two dataset are new ready to be fixed about redundancy, this is removed before splitting for clean these set. Fristly we have merged these two and for remove nearly identical sequences, idenfity representative sequences, build non-redundant tranining sets and filter similar sequences in these. The command used is: 
```
seqs easy-cluster input.fa cluster-results tmp \ --min-seq-id 0.3 -c 0.4 --cov-mode 0 --cluster-mode 1
```
That says we have used a 30% sequence-identity threshold and a 40% length-coverage threshold. 


After that we splitted randomly this set in two further sets each, specifically with 80% in one called training set and the remaining (20%) in other one called test set.
At this point we have two training and two test set, these latter ones will be used in the end. In the two training set we have done the 5-fold CV, that means split the set into 5 subsets, this is done for evaluating the model's performance more reliably.
For preserve the positive/negative ratio in each subset 


In the file data_collection.py there's the entire pipeline, for use it from the terminal the usage is:
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


