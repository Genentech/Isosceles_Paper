#!/usr/bin/env Rscript

suppressMessages({
    library(tidyverse)
    library(glue)
    library(Isosceles)
})

# Set the number of CPUs/threads for the analysis
ncpu <- 1

# Global parameters
genome_fasta_file <- file.path("..", "reference_data", "genome.fasta")
gtf_file <- file.path("..", "reference_data", "simulated_transcript_annotations.gtf")
cell_line_ids <- c("SK-OV-3", "IGROV-1", "OVMANA", "OVKATE", "OVTOKO", "COV362")
sample_ids <- glue("{cell_line_ids}_Rep2")
bam_files <- as.character(glue("bam/{sample_ids}.bam"))
names(bam_files) <- sample_ids

result_dir <- "isosceles_results_bulk"
dir.create(result_dir, recursive = TRUE)

# Preparing transcript set
min_intron_length <- 30
max_intron_length <- 5e6
known_intron_motifs <- c("GT-AG")
rescue_annotated_introns <- TRUE
min_bam_splice_read_count <- 2
min_bam_splice_fraction <- 0.1
bin_size <- 50
transcript_data <- prepare_transcripts(
    gtf_file = gtf_file,
    genome_fasta_file = genome_fasta_file,
    bam_parsed = NULL,
    min_intron_length = min_intron_length,
    max_intron_length = max_intron_length,
    known_intron_motifs = known_intron_motifs,
    rescue_annotated_introns = rescue_annotated_introns,
    known_intron_granges = NULL,
    min_bam_splice_read_count = min_bam_splice_read_count,
    min_bam_splice_fraction = min_bam_splice_fraction,
    bin_size = bin_size
)
saveRDS(transcript_data, file.path(
    result_dir, "transcript_data.rds"
))

# Preparing the TCC SE object
run_mode <- "strict"
min_read_count <- 1
min_relative_expression <- 0
extend_spliced_transcripts <- 100
chunk_size <- 1000000
for (sample_id in sample_ids) {
    bam_file <- bam_files[sample_id]
    se_tcc <- bam_to_tcc(
        bam_files = bam_file,
        transcript_data = transcript_data,
        run_mode = run_mode,
        min_read_count = min_read_count,
        min_relative_expression = min_relative_expression,
        extend_spliced_transcripts = extend_spliced_transcripts,
        chunk_size = chunk_size,
        ncpu = ncpu
    )
    saveRDS(se_tcc, file.path(result_dir, glue("{sample_id}_se_tcc.rds")))
}

# Preparing the transcript SE object
em.maxiter <- 250
em.conv <- 0.01
use_length_normalization <- TRUE
for (sample_id in sample_ids) {
    se_tcc <- readRDS(file.path(result_dir, glue("{sample_id}_se_tcc.rds")))
    se_transcript <- tcc_to_transcript(
        se_tcc = se_tcc,
        em.maxiter = em.maxiter,
        em.conv = em.conv,
        use_length_normalization = use_length_normalization,
        ncpu = ncpu
    )
    saveRDS(se_transcript, file.path(result_dir, glue("{sample_id}_se_transcript.rds")))
}

# Print session information
sessionInfo()
