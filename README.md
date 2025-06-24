# AF3 JSON Generator

This repository contains R scripts for generating AlphaFold3 input JSON files from a list of protein complexes and local FASTA sequences.

## Usage

```bash
Rscript generate_json_AF3.R \
  -f ./human_seqs.fasta \
  -c ./test_complexes.csv \
  -b 30 \
  -o ./output_dir \
  -r ./failed_complexes.csv
```

### Input format
The input file should be a `.csv`, where the **1st** column is the name of the target complex. For example:

```
Complex
P00734-P00533
P00734-O75015
P00734-P00736
P00734-P02745
```
Just so you know, this version supports names with **UNIPROT** ID only, and should be separated by '-'. 


### Protein sequence
The current `human_seqs.fasta` contains the protein sequence of human from [UniProt Reviewed (Swiss-Prot)](https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz). You can also use your own `.fasta` by changing the `-f` parameter.



