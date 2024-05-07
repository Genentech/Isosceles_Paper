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
bam_file <- file.path("wf_single_cell_results", "umi_tools_dedup",
                      "LIB5445493_SAM24404003.bam")
sample_id <- "LIB5445493_SAM24404003"
names(bam_file) <- sample_id
souporcell_file <- file.path("..", "illumina_sc_analysis", "data",
                             "souporcell_clusters.tsv")

result_dir <- "isosceles_results_wf"
dir.create(result_dir, recursive = TRUE)

# Preparing transcript set
min_intron_length <- 30
max_intron_length <- 5e6
known_intron_motifs <- c("GT-AG")
rescue_annotated_introns <- TRUE
min_bam_intron_read_count <- 2
bin_size <- 50
known_intron_granges <- rtracklayer::import(intron_bed_file)
transcript_data <- prepare_transcripts(
    gtf_file = gtf_file, genome_fasta_file = genome_fasta_file,
    bam_parsed = NULL,
    min_intron_length = min_intron_length,
    max_intron_length = max_intron_length,
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
is_single_cell <- TRUE
barcode_tag <- "CB"
chunk_size <- 1000000
se_tcc <- bam_to_tcc(
    bam_files = bam_file, transcript_data = transcript_data,
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
se_gene <- tcc_to_gene(se_tcc = se_tcc)
saveRDS(se_gene, file.path(result_dir, "se_gene.rds"))

# Preparing the transcript SE object
em.maxiter <- 250
em.conv <- 0.01
use_length_normalization <- FALSE
se_transcript <- tcc_to_transcript(
    se_tcc = se_tcc, em.maxiter = em.maxiter, em.conv = em.conv,
    use_length_normalization = use_length_normalization, ncpu = ncpu
)
saveRDS(se_transcript, file.path(result_dir, "se_transcript.rds"))

# Preparing the Souporcell cluster data
souporcell_df <- read.delim(souporcell_file)
souporcell_barcode <- gsub("-1$", "", souporcell_df$barcode)
souporcell_barcode <- paste0(sample_id, ".", souporcell_barcode)
souporcell_cluster <- ifelse(souporcell_df$status == "singlet",
                             souporcell_df$assignment,
                             souporcell_df$status)
souporcell_cluster <- setNames(souporcell_cluster, souporcell_barcode)
cell_labels <- souporcell_cluster[colnames(se_tcc)]
se_tcc$label <- cell_labels
se_gene$label <- cell_labels
se_transcript$label <- cell_labels
saveRDS(se_tcc, file.path(result_dir, "se_tcc.rds"))
saveRDS(se_gene, file.path(result_dir, "se_gene.rds"))
saveRDS(se_transcript, file.path(result_dir, "se_transcript.rds"))

# Preparing the pseudobulk TCC SE object
cell_labels <- se_tcc$label
se_pseudobulk_tcc <- pseudobulk_tcc(se_tcc, cell_labels)
saveRDS(se_pseudobulk_tcc, file.path(result_dir, "se_pseudobulk_tcc.rds"))

# Preparing the pseudobulk gene SE object
se_pseudobulk_gene <- tcc_to_gene(se_tcc = se_pseudobulk_tcc)
saveRDS(se_pseudobulk_gene, file.path(result_dir, "se_pseudobulk_gene.rds"))

# Preparing the pseudobulk transcript SE object
em.maxiter <- 250
em.conv <- 0.01
use_length_normalization <- FALSE
se_pseudobulk_transcript <- tcc_to_transcript(
    se_tcc = se_pseudobulk_tcc, em.maxiter = em.maxiter, em.conv = em.conv,
    use_length_normalization = use_length_normalization, ncpu = ncpu
)
saveRDS(se_pseudobulk_transcript, file.path(result_dir, "se_pseudobulk_transcript.rds"))

# Print session information
sessionInfo()
