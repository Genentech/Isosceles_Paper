#!/usr/bin/env Rscript

suppressMessages({
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
genome_fasta_file <- file.path("reference_data", "sirvome.fasta")
annotations_ids <- c("insufficient", "over")
gtf_files <- as.character(glue("reference_data/sirvome_{annotations_ids}_annotations.gtf"))
names(gtf_files) <- annotations_ids
sample_id <- "cdna"
bam_file <- as.character(glue("bam/sirvome_{sample_id}.bam"))
names(bam_file) <- sample_id

# Run bambu
for (annotations_id in annotations_ids) {
    gtf_file <- gtf_files[annotations_id]
    txdb <- suppressWarnings(makeTxDbFromGFF(gtf_file))
    annotations <- prepareAnnotations(txdb)
    se_bambu <- bambu(reads = bam_file, annotations = annotations,
                      genome = genome_fasta_file, discovery = TRUE,
                      ncore = ncpu, yieldSize = 1e6, lowMemory = TRUE)
    saveRDS(se_bambu, file.path(bambu_dir, glue("{sample_id}_{annotations_id}.rds")))
    writeBambuOutput(se_bambu, path = bambu_dir, prefix = glue("{sample_id}_{annotations_id}_"))
}

# Print session information
sessionInfo()
