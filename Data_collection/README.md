# Data collection

`data_collection.py` downloads the positive and the negative set from the UniProtKB REST API and writes each of them as a TSV table and a FASTA file.

```bash
python3 data_collection.py -pos <output prefix>
python3 data_collection.py -neg <output prefix>
python3 data_collection.py -all <output prefix negative> <output prefix positive>
```

Each prefix produces `<prefix>.tsv` and `<prefix>.fasta`. The files in this folder were created, from the repository root, with:

```bash
python3 Data_collection/data_collection.py -all Data_collection/negative/negative_set Data_collection/positive_set
```

Entries are downloaded in batches of 500, following the pagination links returned by UniProt; requests that fail with HTTP 500, 502, 503 or 504 are retried up to 5 times.

## Positive set

UniProtKB query:

```
(ft_signal_exp:*) AND (reviewed:true) AND (fragment:false) AND (taxonomy_id:2759) AND (length:[40 TO *]) AND (existence:1)
```

that is, reviewed, non-fragment eukaryotic proteins of at least 40 residues, with evidence at protein level and an experimentally verified signal peptide. The query returns 2974 entries. The script then keeps only the proteins whose signal peptide starts at residue 1, is at least 14 residues long and is not annotated as "Not Cleaved": **2961 proteins** (`positive_set.tsv`, `positive_set.fasta`). The cleavage position is the last residue of that same signal peptide.

## Negative set

UniProtKB query:

```
(taxonomy_id:2759) AND (existence:1) AND (reviewed:true) AND (fragment:false) AND (length:[40 TO *])
AND ((cc_scl_term_exp:SL-0091) OR (cc_scl_term_exp:SL-0191) OR (cc_scl_term_exp:SL-0173) OR (cc_scl_term_exp:SL-0039) OR (cc_scl_term_exp:SL-0209) OR (cc_scl_term_exp:SL-0204))
NOT (ft_signal:*)
NOT ((cc_scl_term:SL-0095) OR (cc_scl_term:SL-0132) OR (cc_scl_term:SL-0112) OR (cc_scl_term:SL-0243))
```

that is, the same constraints as the positives, no signal peptide at any evidence level, and an experimentally verified localization in cytosol, nucleus, mitochondrion, plastid, peroxisome or cell membrane. Proteins also annotated, at any evidence level, in the endoplasmic reticulum, Golgi apparatus, extracellular space or as secreted are excluded: this removes 1380 of the 20974 entries returned by the first version of the query, leaving **19594 proteins** (`negative/negative_set.tsv`, `negative/negative_set.fasta`). No further filter is applied.

For each negative we also record whether a helical transmembrane segment starts before residue 90 (`early_trans`), because a hydrophobic transmembrane helix near the N-terminus can be mistaken for a signal peptide.

## Output format

TSV without header, one protein per row:

| Column | Content |
|---|---|
| 1 | class: 1 = positive, 0 = negative |
| 2 | UniProt accession |
| 3 | organism name |
| 4 | kingdom: Metazoa, Fungi, Plants (Viridiplantae) or Other |
| 5 | protein length |
| 6 | cleavage position, i.e. the last residue of the signal peptide (0 for negatives) |
| 7 | `early_trans`: helical transmembrane segment starting before residue 90 (always False for positives) |

The FASTA file has the UniProt accession as header and each sequence on a single line.

## Main functions

| Function | Role |
|---|---|
| `get_next_link`, `get_batch` | follow the UniProt pagination and yield one batch of results at a time |
| `valid_signal` | returns the signal peptide that satisfies the positive criteria, or None |
| `filter_function_positive` | True if `valid_signal` finds a signal peptide |
| `extract_fields_positive`, `extract_fields_negative` | extract the TSV fields and the sequence from a UniProt JSON entry |
| `get_dataset(url, extract_function, output_prefix, [filter_function])` | downloads all the batches, applies the filter if given and writes the TSV and FASTA files |

## Older scripts

`positive/positive_filter.py`, `negative/negative_dataset.py` and their notebooks are the first versions of the download, kept for reference. They write a 5-column TSV without the class column, which `remove_redundancy.py` and `datasplitting.py` cannot read: use `data_collection.py`.
