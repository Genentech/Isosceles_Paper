#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(GenomicFeatures)
library(rtracklayer)
library(Biostrings)
library(BSgenome)

# Global parameters
genome_fasta <- file.path("reference_data", "sequin.fasta")
gtf_file <- file.path("reference_data", "sequin_annotations.gtf")

# Read the reference genome sequence
genome_seq <- readDNAStringSet(genome_fasta, format = "fasta")

# Read genome annotations
txdb <- suppressWarnings(makeTxDbFromGFF(gtf_file))

# Prepare the reference transcriptome sequence
transcriptome_seq <- extractTranscriptSeqs(genome_seq, txdb, use.names = TRUE)
writeXStringSet(transcriptome_seq,
                file.path("reference_data", "transcriptome.fasta"),
                format = "fasta")

# Prepare intron BED files
export(
    unique(tidyIntrons(txdb)),
    file.path("reference_data", "known_introns.bed"),
    format = "bed"
)

# Print session information
sessionInfo()
