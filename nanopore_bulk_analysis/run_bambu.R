#!/usr/bin/env Rscript

suppressMessages({
  library(devtools)
  library(tidyverse)
  library(glue)
  library(GenomicFeatures)
  library(bambu)
})

# Set the number of CPUs/threads for the analysis
ncpu <- 1

# Global parameters
bambu_dir <- "bambu_results"
dir.create(bambu_dir, recursive = TRUE)
genome_fasta_file <- file.path("..", "reference_data", "genome.fasta")
sample_ids <- c("LIB5432315_SAM24385458", "LIB5432316_SAM24385459",
                "LIB5427896_SAM24376275", "LIB5427897_SAM24376276")
bam_files <- as.character(glue("bam/{sample_ids}.bam"))
names(bam_files) <- sample_ids

# Run bambu
gtf_file <- file.path("..", "reference_data", "ensembl_90.gtf")
txdb <- suppressWarnings(makeTxDbFromGFF(gtf_file))
annotations <- prepareAnnotations(txdb)
for (sample_id in sample_ids) {
    bam_file <- bam_files[sample_id]
    se_bambu <- bambu(reads = bam_file, annotations = annotations,
                      genome = genome_fasta_file, discovery = FALSE,
                      ncore = ncpu, yieldSize = 1e6, lowMemory = TRUE)
    saveRDS(se_bambu, file.path(bambu_dir, glue("{sample_id}.rds")))
    writeBambuOutput(se_bambu, path = bambu_dir, prefix = glue("{sample_id}_"))
}

# Print session information
sessionInfo()
