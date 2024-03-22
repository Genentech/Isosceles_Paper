#!/usr/bin/env Rscript

suppressMessages({
    library(tidyverse)
    library(glue)
    library(Isosceles)
})

# Set the number of CPUs/threads for the analysis
ncpu <- 1

# Global parameters
genome_fasta_file <- file.path("reference_data", "sirvome.fasta")
sample_id <- "cdna"
bam_file <- as.character(glue("bam/sirvome_{sample_id}.bam"))
names(bam_file) <- sample_id
all_annotations_ids <- c("correct", "insufficient", "over")
annotations_ids <- c("insufficient", "over")
gtf_files <- as.character(glue("reference_data/sirvome_{all_annotations_ids}_annotations.gtf"))
names(gtf_files) <- all_annotations_ids

result_dir <- "isosceles_results"
dir.create(result_dir, recursive = TRUE)

# Prepare reference annotation data
dir.create(file.path(result_dir, "reference_data"), recursive = TRUE)
bin_size <- 50
for (annotations_id in all_annotations_ids) {
    gtf_file <- gtf_files[annotations_id]
    transcript_data <- prepare_transcripts(
        gtf_file = gtf_file,
        genome_fasta_file = genome_fasta_file,
        bam_parsed = NULL,
        bin_size = bin_size
    )
    saveRDS(transcript_data, file.path(
        result_dir, "reference_data",
        glue("annotations_{annotations_id}.rds")))
}

# BAM parsing
dir.create(file.path(result_dir, "bam_parsed"), recursive = TRUE)
bam_parsed <- bam_to_read_structures(bam_file)
saveRDS(bam_parsed, file.path(result_dir, "bam_parsed",
                              glue("{sample_id}.rds")))

# Preparing transcript set
dir.create(file.path(result_dir, "transcript_data"), recursive = TRUE)
min_intron_length <- 30
max_intron_length <- 5e6
known_intron_motifs <- c("GT-AG")
rescue_annotated_introns <- TRUE
min_bam_splice_read_count <- 2
min_bam_splice_fraction <- 0.1
bin_size <- 50
for (annotations_id in annotations_ids) {
    gtf_file <- gtf_files[annotations_id]
    transcript_data <- prepare_transcripts(
        gtf_file = gtf_file,
        genome_fasta_file = genome_fasta_file,
        bam_parsed = bam_parsed,
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
        result_dir, "transcript_data",
        glue("{sample_id}_{annotations_id}.rds")
    ))
}

# Preparing the TCC SE objects
run_mode <- "de_novo_loose"
min_read_count <- 50
min_relative_expression <- 0
extend_spliced_transcripts <- 100
chunk_size <- 1000000
for (annotations_id in annotations_ids) {
    transcript_data <- readRDS(file.path(
        result_dir, "transcript_data",
        glue("{sample_id}_{annotations_id}.rds")
    ))
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
    saveRDS(se_tcc, file.path(
        result_dir, glue("{sample_id}_{annotations_id}_se_tcc.rds")
    ))
}

# Preparing the transcript SE objects
em.maxiter <- 250
em.conv <- 0.01
use_length_normalization <- TRUE
for (annotations_id in annotations_ids) {
    se_tcc <- readRDS(file.path(
        result_dir, glue("{sample_id}_{annotations_id}_se_tcc.rds")
    ))
    se_transcript <- tcc_to_transcript(
        se_tcc = se_tcc,
        em.maxiter = em.maxiter,
        em.conv = em.conv,
        use_length_normalization = use_length_normalization,
        ncpu = ncpu
    )
    saveRDS(se_transcript, file.path(
        result_dir, glue("{sample_id}_{annotations_id}_se_transcript.rds")
    ))
}

# Print session information
sessionInfo()
