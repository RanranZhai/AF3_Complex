#!/usr/bin/env Rscript

# 自动安装并加载包
ensure_packages <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, dependencies = TRUE)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}
ensure_packages(c('argparse', 'Biostrings', 'stringr'))


## create JSON entity for a single protein sequence
entity.ot <- function(seq) {
	(paste0(
      '\n{
        \"proteinChain\": {
          \"sequence\": \"', seq, '\",
          \"count\": 1
        }
      },\n'))
}


## Create JSON file for list of complex in batch mode
crt.json.batch <- function(job.name, seqs, counts) {
	ot <- ''
	for (seq in seqs) {
		ot <- paste0(ot, entity.ot(seq))
	}
	ot <- str_sub(ot, end = -3)
	paste0('
  {
    \"name\": \"', job.name, '\",
    \"modelSeeds\": [],
    \"sequences\": [',
      ot,
    ']
  },'
	)
}

crt.json.batch.last <- function(job.name, seqs, counts) {
	ot <- ''
	for (seq in seqs) {
		ot <- paste0(ot, entity.ot(seq))
	}
	ot <- str_sub(ot, end = -3)
	paste0('
  {
    \"name\": \"', job.name, '\",
    \"modelSeeds\": [],
    \"sequences\": [',
      ot,
    ']
  }'
	)
}


library(argparse)
library(Biostrings)
library(stringr)

# arguments
parser <- ArgumentParser(description = 'Generate AF3 JSON files from protein complexes and local fasta sequences')
parser$add_argument('-f', '--fasta', default = './human_seqs.fasta', help = 'Path to local FASTA file containing sequences')
parser$add_argument('-c', '--complex', required = TRUE, help = 'CSV file listing complexes ("-" separated)')
parser$add_argument('-b', '--batchsize',  default = 30, help = 'number of complexes per JSON file (default: 30)')
#parser$add_argument('--id_type', default = 'uniprot', help = 'Use "uniprot" (default) or "gene" as sequence names')
parser$add_argument('-o', '--outdir', default = 'AF3_JSON', help = 'Output directory for JSON files')
parser$add_argument('-r', '--report', default = 'failed_complexes.csv', help = 'CSV report filename for failed complexes')
args <- parser$parse_args()


# read fasta file
all_seqs <- readAAStringSet(args$fasta)
#all_seqs <- readAAStringSet('human_seqs.fasta')
names(all_seqs) <- toupper(sapply(strsplit(names(all_seqs), "\\|"), `[`, 2))


## if (id_type == 'uniprot') {
##   names(all_seqs) <- toupper(sapply(strsplit(names(all_seqs), "\\|"), `[`, 2))
## } else if (id_type == 'gene') {
##   names(all_seqs) <- sapply(stringr::str_match(names(all_seqs), 'GN=([^ ]+)'), `[`, 2)
## } else {
##   stop("Invalid id_type. Choose from: 'uniprot' or 'gene'")
## }
## head(names(all_seqs))


# read complex list
complex_df <- read.csv(args$complex, stringsAsFactors = FALSE)
cmps <- complex_df[, 1]
uniprot_ids <- unique(unlist(strsplit(cmps, "-")))


missing_cmps <- missing_pros <- c()
for (cmp in cmps) {
    ps <- unlist(strsplit(cmp, '-'))
    if (sum(ps %in% names(all_seqs)) != length(ps)) {
        missing_cmps <- c(missing_cmps, cmp)
        #missing_pros <- c(missing_pros, ps[!ps %in% names(all_seqs)])
    }
}
missing_cmps <- unique(missing_cmps)

write.csv(data.frame(MissingComplexes = missing_cmps), file = 'failed_complexes.csv', row.names = FALSE)

if (length(missing_cmps) > 0) {
    warning('The following complexes failed due to missing in the FASTA file, Please check the ', args$report)
    write.csv(data.frame(MissingComplexes = missing_cmps), file = args$report, row.names = FALSE)
}


cmps <- cmps[!cmps %in% missing_cmps]  # remove missing complexes from the list
message('Creating JSON files for ', length(cmps), ' complexes...')


## output json files in batch mode
num <- as.numeric(args$batchsize)
N <- ceiling(length(cmps)/num)

if(!dir.exists(args$outdir)) {
    dir.create(args$outdir, recursive = TRUE)
}

for (i in 1:N) {
    a <- (i-1)*num+1
    b <- min((i*num), length(cmps))

    json_content <- '[\n'
    for (k in a:b) {
        cmp <- cmps[k]
        gs <- unlist(strsplit(cmp, '-'))
	    seq <- list()
	    for (g in gs) {
			seq[[g]] <- as.character(all_seqs[g])
	    }
        if (k == b) {
            json_content <- paste0(json_content, crt.json.batch.last(cmp, seq, rep(1, times = length(seq))))
        } else {
            json_content <- paste0(json_content, crt.json.batch(cmp, seq, rep(1, times = length(seq))))
        }
    }
    json_content <- paste0(json_content, ']\n')
    writeLines(json_content, paste0(args$outdir, '/job_batch', i, '.json'))
}

message("✅ JSON files generated successfully in directory: ", args$outdir)


# Running the script:
# Rscript generate_json_AF3.R -f human_seqs.fasta -c test_complexes.csv -b 30 -o AF3_JSON -r failed_complexes.csv
