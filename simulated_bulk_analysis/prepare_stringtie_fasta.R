#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(rtracklayer)
library(GenomicFeatures)
library(Biostrings)
library(BSgenome)

# Global parameters
data_dir <- "../reference_data"
stringtie_dir <- "stringtie_results"
genome_fasta <- file.path(data_dir, "genome.fasta")
perc_downs <- c(10, 20, 30)

# Read the reference genome sequence
genome_seq <- readDNAStringSet(genome_fasta, format = "fasta")
names(genome_seq) <- sapply(strsplit(names(genome_seq), "\\s+"), "[", 1)

# Extract transcript sequences from StringTie GTF files
for (perc_down in perc_downs) {
    gtf_file <- file.path(stringtie_dir, glue("truncated_bulk_rnaseq_denovo_{perc_down}.gtf"))
    txdb <- suppressWarnings(makeTxDbFromGFF(gtf_file))
    transcriptome_seq <- extractTranscriptSeqs(genome_seq, txdb, use.names = TRUE)
    writeXStringSet(transcriptome_seq,
                    file.path(stringtie_dir, glue("truncated_bulk_rnaseq_denovo_{perc_down}.fasta")),
                    format = "fasta")
}

# Print session information
sessionInfo()
