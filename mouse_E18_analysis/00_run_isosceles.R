#!/usr/bin/env Rscript

suppressMessages({
    library(tidyverse)
    library(glue)
    library(Isosceles)
})

# Set the number of CPUs/threads for the analysis
ncpu <- 1

# Global parameters
genome_fasta_file <- file.path("data", "genome.fasta")
gtf_file <- file.path("data", "gencode_m25.gtf")
intron_bed_file <- file.path("data", "known_introns.bed")
sample_ids <- c("190_cells", "951_cells")
bam_files <- glue("data/{sample_ids}.bam")
names(bam_files) <- sample_ids

result_dir <- "00_run_isosceles"
dir.create(result_dir, recursive = TRUE)

# BAM parsing
chunk_size <- 1000000
bam_parsed <- bam_to_read_structures(
    bam_files = bam_files, chunk_size = chunk_size, ncpu = ncpu
)
saveRDS(bam_parsed, file.path(result_dir, "bam_parsed.rds"))

# Preparing transcript set
min_intron_length <- 30
max_intron_length <- 5e6
known_intron_motifs <- c("GT-AG")
rescue_annotated_introns <- FALSE
min_bam_intron_read_count <- 2
bin_size <- 50
known_intron_granges <- rtracklayer::import(intron_bed_file)
bam_parsed <- readRDS(file.path(result_dir, "bam_parsed.rds"))
transcript_data <- prepare_transcripts(
    gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
    bam_parsed = bam_parsed,
    min_intron_length = min_intron_length,
    max_intron_length = max_intron_length,
    known_intron_motifs = known_intron_motifs,
    rescue_annotated_introns = rescue_annotated_introns,
    known_intron_granges = known_intron_granges,
    bin_size = bin_size
)
saveRDS(transcript_data, file.path(result_dir, "transcript_data.rds"))

# Preparing the TCC SE object
run_mode <- "de_novo_loose"
min_read_count <- 1
min_relative_expression <- 0
extend_spliced_transcripts <- 100
is_single_cell <- TRUE
barcode_tag <- "BC"
chunk_size <- 1000000
transcript_data <- readRDS(file.path(result_dir, "transcript_data.rds"))
se_tcc <- bam_to_tcc(
    bam_files = bam_files, transcript_data = transcript_data,
    run_mode = run_mode,
    min_read_count = min_read_count,
    min_relative_expression = min_relative_expression,
    extend_spliced_transcripts = extend_spliced_transcripts,
    is_single_cell = is_single_cell,
    barcode_tag = barcode_tag,
    chunk_size = chunk_size, ncpu = ncpu
)
saveRDS(se_tcc, file.path(result_dir, "se_tcc.rds"))

# Preparing the gene SE object
se_tcc <- readRDS(file.path(result_dir, "se_tcc.rds"))
se_gene <- tcc_to_gene(se_tcc = se_tcc)
saveRDS(se_gene, file.path(result_dir, "se_gene.rds"))

# Preparing the transcript SE object
em.maxiter <- 250
em.conv <- 0.01
use_length_normalization <- FALSE
se_tcc <- readRDS(file.path(result_dir, "se_tcc.rds"))
se_transcript <- tcc_to_transcript(
    se_tcc = se_tcc, em.maxiter = em.maxiter, em.conv = em.conv,
    use_length_normalization = use_length_normalization, ncpu = ncpu
)
saveRDS(se_transcript, file.path(result_dir, "se_transcript.rds"))

# Preparing the PSI SE object
se_transcript <- readRDS(file.path(result_dir, "se_transcript.rds"))
se_psi <- transcript_to_psi(se_transcript, ncpu = ncpu)
saveRDS(se_psi, file.path(result_dir, "se_psi.rds"))

# Print session information
sessionInfo()
