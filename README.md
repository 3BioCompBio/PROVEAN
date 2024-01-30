# PROVEAN (Protein Variation Effect Analyzer) Python implementation

This `provean.py` file provides a Python implementation* of the famous [PROVEAN method](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0046688), which predicts the functional effect of protein variants, including amino acid substitutions and indels. 

## Dependencies

Requires `biopython>=1.80` and `joblib>=1.3.2` to be installed, as well as CD-HIT.

Unlike the implementation distributed by the original authors, this implementation does not require you to install a sequence database such as NR and run BLAST on your machine, and instead uses the NCBI REST API to run BLASTP on their servers to obtain homologous sequences. Also, if you have your own pipeline or a pre-computed fasta file with homologous sequences, you can provide it through the `-homologous_sequences_fasta_file` argument and avoid running BLASTP completely.  

## Usage

The command line interface takes as input a FASTA with your protein of interest (1 chain only) and an output directory path. For example, to run PROVEAN on the human protein [P53](https://www.uniprot.org/uniprotkb/P04637/entry#sequences), the command would look something like

```
python provean.py -n_cores=4 /home/user/Downloads/P04637.fasta /home/user/output_dir
```

The file with the predicted PROVEAN scores will be written in the output directory as "PROVEAN_scores.csv", and intermediate files such as the list of homologous sequences returned by the BLASTP or the clusters file from CD-HIT can also be obtained by setting `-save_intermediate_files=True`. To see all the available arguments and parameters, run `python provean.py -h`.

Regarding the interpretation of the scores predicted by PROVEAN, the original authors recommend to use a threshold of -2.5, where values below this threshold are predicted as deleterious, and those above are predicted to be neutral.

***Note: when predicting the PROVEAN score for amino acid substitutions, our implementation does not realign the query sequence variants against the selected homologous sequences, and instead only computes the alignment between the wild type query sequence and the homologous sequences once and then uses the BLOSUM62 matrix to rescore the alignment for the given amino acid substitution. The original authors perform a partial realignment of each sequence variant against the homologous sequences, but we have not re-implemented this, and performing a full realingment for each variant would be very slow. As a result, our implementation tends to predict slightly more negative scores than the original implementation. For indels, we do perform a full realignment, and therefore these scores do not have this bias.**




