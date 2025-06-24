# AF3 JSON Generator

This repository contains R scripts for generating AlphaFold3 input JSON files from a list of protein complexes and local FASTA sequences.

## Usage

```bash
Rscript generate_json_AF3.R \
  -f ./human_seqs.fasta \
  -c ./test_complexes.csv \
  -o ./output_dir \
  -r ./failed_complexes.csv
