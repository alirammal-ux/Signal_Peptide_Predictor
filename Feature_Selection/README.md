## Feature Extraction and Selection:


Features are extracted using chemico-phisical scales from ProtScale and AAindex. The scales are applied through sliding window to the first 40 residues of each protein sequence.

| scale | window size | motivation |
| :--- | :---: |:--- |
| Hydropathy index (Kyte-Doolittle) | 5 | h-region in signal peptide is hydrophobic (6,15 residues) |
| Positive charge (Fauchere) | [3,5] | n-region in signal peptide is positivelly charged while the c-region of the cleavage site contains uncharged residues |
| Negative charge (Faucher) | [3.5] |
| Normalized frequency of N-terminal helix (Chou-Fasman) | [7,11] | h-region in signal peptide contains alpha helix (6,15 residues) |
| Residue volume (Bigelo)| [3,5] | c-region in signal peptide contains small residues |
| Signal sequence helical potential (Argos) | [7,11] | this scale quantifies the likelihood of an amino acid residue to adopt an alpha-helix conformation within the context of a signal peptide |
| Average flexibility indices (Bhaskaran-Ponnuswamy) | [7,11] | 
| Transmembrane tendency | [7,11] | to discriminate between sequences that start with a hydrophobic region due to signal peptides and sequences for transmembran protein |

From the hydrobhobicity profile of each sequence, we extract the following features:
- hydrophobicity score of the firs maximum peak
- position in the sequence of the first maximum peak
- length (= number of positions) of the region around the peak with hydrophobicity score upper a threshold
- drop in hydrophobicity = hydrophobicity score of the peak minus the hidrophobicity score of the 40° residue

Other 20 features are computed as the frequency of each residue in the first n residues in the protein sequence, where n is the median SP length in the training set (recomputed on the training folds of each CV split).


