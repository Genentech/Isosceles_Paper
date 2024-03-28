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
genome_fasta_file <- file.path("..", "reference_data", "genome.fasta")
sample_id <- "truncated_bulk_rnaseq"
bam_file <- as.character(glue("bam/{sample_id}.bam"))
names(bam_file) <- sample_id
perc_downs <- c(10, 20, 30)

# Run bambu (transcript quantification)
gtf_file <- file.path("..", "reference_data", "benchmark_transcript_annotations.gtf")
txdb <- suppressWarnings(makeTxDbFromGFF(gtf_file))
annotations <- prepareAnnotations(txdb)

se_bambu <- bambu(reads = bam_file, annotations = annotations,
                  genome = genome_fasta_file, discovery = FALSE,
                  ncore = ncpu, yieldSize = 1e6, lowMemory = TRUE)
saveRDS(se_bambu, file.path(bambu_dir, glue("{sample_id}_quant.rds")))
writeBambuOutput(se_bambu, path = bambu_dir, prefix = glue("{sample_id}_quant_"))

# Run bambu (de novo detection)
for (perc_down in perc_downs) {
    gtf_file <- file.path("..", "reference_data", glue("benchmark_downsampled_{perc_down}.gtf"))
    txdb <- suppressWarnings(makeTxDbFromGFF(gtf_file))
    annotations <- prepareAnnotations(txdb)

    se_bambu <- bambu(reads = bam_file, annotations = annotations,
                      genome = genome_fasta_file, discovery = TRUE,
                      ncore = ncpu, yieldSize = 1e6, lowMemory = TRUE)
    saveRDS(se_bambu, file.path(bambu_dir, glue("{sample_id}_denovo_{perc_down}.rds")))
    writeBambuOutput(se_bambu, path = bambu_dir, prefix = glue("{sample_id}_denovo_{perc_down}_"))
}

# Run bambu (de novo detection using StringTie)
for (perc_down in perc_downs) {
    gtf_file <- file.path("stringtie_results", glue("{sample_id}_denovo_{perc_down}.gtf"))
    txdb <- suppressWarnings(makeTxDbFromGFF(gtf_file))
    annotations <- prepareAnnotations(txdb)

    se_bambu <- bambu(reads = bam_file, annotations = annotations,
                      genome = genome_fasta_file, discovery = TRUE,
                      ncore = ncpu, yieldSize = 1e6, lowMemory = TRUE)
    saveRDS(se_bambu, file.path(bambu_dir, glue("{sample_id}_stringtie_{perc_down}.rds")))
    writeBambuOutput(se_bambu, path = bambu_dir, prefix = glue("{sample_id}_stringtie_{perc_down}_"))
}

# Print session information
sessionInfo()
