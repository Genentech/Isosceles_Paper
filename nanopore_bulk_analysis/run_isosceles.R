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
gtf_file <- file.path("..", "reference_data", "ensembl_90.gtf")
intron_bed_file <- file.path("..", "reference_data", "known_introns.bed")
sample_ids <- c("LIB5432309_SAM24385452", "LIB5432310_SAM24385453",
                "LIB5432311_SAM24385454", "LIB5432312_SAM24385455",
                "LIB5432313_SAM24385456", "LIB5432314_SAM24385457",
                "LIB5432315_SAM24385458", "LIB5432316_SAM24385459")
bam_files <- as.character(glue("bam/{sample_ids}.bam"))
names(bam_files) <- sample_ids
sample_ids_2 <- c("LIB5432315_SAM24385458", "LIB5432316_SAM24385459",
                  "LIB5427896_SAM24376275", "LIB5427897_SAM24376276")
bam_files_2 <- as.character(glue("bam/{sample_ids_2}.bam"))
names(bam_files_2) <- sample_ids_2

result_dir <- "isosceles_results"
dir.create(result_dir, recursive = TRUE)

# Preparing transcript set
min_intron_length <- 30
known_intron_motifs <- c("GT-AG")
rescue_annotated_introns <- FALSE
min_bam_intron_read_count <- 2
bin_size <- 50
known_intron_granges <- rtracklayer::import(intron_bed_file)
transcript_data <- prepare_transcripts(
    gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
    bam_parsed = NULL,
    min_intron_length = min_intron_length,
    known_intron_motifs = known_intron_motifs,
    rescue_annotated_introns = rescue_annotated_introns,
    known_intron_granges = known_intron_granges,
    bin_size = bin_size
)
saveRDS(transcript_data, file.path(result_dir, "transcript_data.rds"))

# Preparing the TCC SE object
run_mode <- "strict"
min_read_count <- 1
min_relative_expression <- 0
extend_spliced_transcripts <- 100
chunk_size <- 1000000
transcript_data <- readRDS(file.path(result_dir, "transcript_data.rds"))
se_tcc <- bam_to_tcc(
    bam_files = bam_files, transcript_data = transcript_data,
    run_mode = run_mode,
    min_read_count = min_read_count,
    min_relative_expression = min_relative_expression,
    extend_spliced_transcripts = extend_spliced_transcripts,
    chunk_size = chunk_size, ncpu = ncpu
)
saveRDS(se_tcc, file.path(result_dir, "se_tcc.rds"))
for (sample_id in sample_ids_2) {
    bam_file <- bam_files_2[sample_id]
    se_tcc <- bam_to_tcc(
        bam_files = bam_file, transcript_data = transcript_data,
        run_mode = run_mode,
        min_read_count = min_read_count,
        min_relative_expression = min_relative_expression,
        extend_spliced_transcripts = extend_spliced_transcripts,
        chunk_size = chunk_size, ncpu = ncpu
    )
    saveRDS(se_tcc, file.path(result_dir, glue("{sample_id}_se_tcc.rds")))
}

# Preparing the transcript SE object
em.maxiter <- 250
em.conv <- 0.01
use_length_normalization <- TRUE
se_tcc <- readRDS(file.path(result_dir, "se_tcc.rds"))
se_transcript <- tcc_to_transcript(
    se_tcc = se_tcc, em.maxiter = em.maxiter, em.conv = em.conv,
    use_length_normalization = use_length_normalization, ncpu = ncpu
)
saveRDS(se_transcript, file.path(result_dir, "se_transcript.rds"))
for (sample_id in sample_ids_2) {
    se_tcc <- readRDS(file.path(result_dir, glue("{sample_id}_se_tcc.rds")))
    se_transcript <- tcc_to_transcript(
        se_tcc = se_tcc, em.maxiter = em.maxiter, em.conv = em.conv,
        use_length_normalization = use_length_normalization, ncpu = ncpu
    )
    saveRDS(se_transcript, file.path(result_dir, glue("{sample_id}_se_transcript.rds")))
}

# Preparing the gene SE object
se_tcc <- readRDS(file.path(result_dir, "se_tcc.rds"))
se_gene <- tcc_to_gene(se_tcc = se_tcc)
saveRDS(se_gene, file.path(result_dir, "se_gene.rds"))
for (sample_id in sample_ids_2) {
    se_tcc <- readRDS(file.path(result_dir, glue("{sample_id}_se_tcc.rds")))
    se_gene <- tcc_to_gene(se_tcc = se_tcc)
    saveRDS(se_gene, file.path(result_dir, glue("{sample_id}_se_gene.rds")))
}

# Print session information
sessionInfo()
