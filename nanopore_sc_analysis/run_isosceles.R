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
bam_file <- file.path("sicelore", "molecules_GE_tags.bam")
sample_id <- "LIB5445493_SAM24404003"
names(bam_file) <- sample_id
souporcell_file <- file.path("..", "illumina_sc_analysis", "data",
                             "souporcell_clusters.tsv")
k_values <- c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 800)

result_dir <- "isosceles_results"
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
barcode_tag <- "BC"
chunk_size <- 1000000
transcript_data <- readRDS(file.path(result_dir, "transcript_data.rds"))
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

# Preparing the Souporcell cluster data
souporcell_df <- read.delim(souporcell_file)
souporcell_barcode <- gsub("-1$", "", souporcell_df$barcode)
souporcell_barcode <- paste0(sample_id, ".", souporcell_barcode)
souporcell_cluster <- ifelse(souporcell_df$status == "singlet",
                             souporcell_df$assignment,
                             souporcell_df$status)
souporcell_cluster <- setNames(souporcell_cluster, souporcell_barcode)
se_tcc <- readRDS(file.path(result_dir, "se_tcc.rds"))
se_gene <- readRDS(file.path(result_dir, "se_gene.rds"))
se_transcript <- readRDS(file.path(result_dir, "se_transcript.rds"))
cell_labels <- souporcell_cluster[colnames(se_tcc)]
se_tcc$label <- cell_labels
se_gene$label <- cell_labels
se_transcript$label <- cell_labels
saveRDS(se_tcc, file.path(result_dir, "se_tcc.rds"))
saveRDS(se_gene, file.path(result_dir, "se_gene.rds"))
saveRDS(se_transcript, file.path(result_dir, "se_transcript.rds"))

# Preparing the pseudobulk TCC SE object
se_tcc <- readRDS(file.path(result_dir, "se_tcc.rds"))
cell_labels <- se_tcc$label
se_pseudobulk_tcc <- pseudobulk_tcc(se_tcc, cell_labels)
saveRDS(se_pseudobulk_tcc, file.path(result_dir, "se_pseudobulk_tcc.rds"))

# Preparing the pseudobulk gene SE object
se_pseudobulk_tcc <- readRDS(file.path(result_dir, "se_pseudobulk_tcc.rds"))
se_pseudobulk_gene <- tcc_to_gene(se_tcc = se_pseudobulk_tcc)
saveRDS(se_pseudobulk_gene, file.path(result_dir, "se_pseudobulk_gene.rds"))

# Preparing the pseudobulk transcript SE object
em.maxiter <- 250
em.conv <- 0.01
use_length_normalization <- FALSE
se_pseudobulk_tcc <- readRDS(file.path(result_dir, "se_pseudobulk_tcc.rds"))
se_pseudobulk_transcript <- tcc_to_transcript(
    se_tcc = se_pseudobulk_tcc, em.maxiter = em.maxiter, em.conv = em.conv,
    use_length_normalization = use_length_normalization, ncpu = ncpu
)
saveRDS(se_pseudobulk_transcript, file.path(result_dir, "se_pseudobulk_transcript.rds"))

# Preparing the IGROV-1 TCC SE object
se_tcc <- readRDS(file.path(result_dir, "se_tcc.rds"))
cell_labels <- se_tcc$label
cell_selector <- cell_labels == "1"
se_igrov_tcc <- se_tcc[, cell_selector]
cell_order <- order(colSums(assay(se_igrov_tcc, "counts")), decreasing = TRUE)
se_igrov_tcc <- se_igrov_tcc[, cell_order]
saveRDS(se_igrov_tcc, file.path(result_dir, "se_igrov_tcc.rds"))

# Preparing the IGROV-1 gene SE object
se_igrov_tcc <- readRDS(file.path(result_dir, "se_igrov_tcc.rds"))
se_igrov_gene <- tcc_to_gene(se_tcc = se_igrov_tcc)
saveRDS(se_igrov_gene, file.path(result_dir, "se_igrov_gene.rds"))

# Preparing the IGROV-1 transcript SE object
se_igrov_tcc <- readRDS(file.path(result_dir, "se_igrov_tcc.rds"))
em.maxiter <- 250
em.conv <- 0.01
use_length_normalization <- FALSE
se_igrov_transcript <- tcc_to_transcript(
    se_tcc = se_igrov_tcc, em.maxiter = em.maxiter, em.conv = em.conv,
    use_length_normalization = use_length_normalization, ncpu = ncpu
)
saveRDS(se_igrov_transcript, file.path(result_dir, "se_igrov_transcript.rds"))

# Preparing the IGROV-1 pseudobulk TCC SE object
se_igrov_tcc <- readRDS(file.path(result_dir, "se_igrov_tcc.rds"))
cell_labels <- rep(paste0(sample_id, "_IGROV"), ncol(se_igrov_tcc))
se_igrov_pseudobulk_tcc <- pseudobulk_tcc(se_igrov_tcc, cell_labels)
saveRDS(se_igrov_pseudobulk_tcc, file.path(result_dir, "se_igrov_pseudobulk_tcc.rds"))

# Preparing the IGROV-1 pseudobulk gene SE object
se_igrov_pseudobulk_tcc <- readRDS(file.path(result_dir, "se_igrov_pseudobulk_tcc.rds"))
se_igrov_pseudobulk_gene <- tcc_to_gene(se_tcc = se_igrov_pseudobulk_tcc)
saveRDS(se_igrov_pseudobulk_gene, file.path(result_dir, "se_igrov_pseudobulk_gene.rds"))

# Preparing the IGROV-1 pseudobulk transcript SE object
se_igrov_pseudobulk_tcc <- readRDS(file.path(result_dir, "se_igrov_pseudobulk_tcc.rds"))
em.maxiter <- 250
em.conv <- 0.01
use_length_normalization <- FALSE
se_igrov_pseudobulk_transcript <- tcc_to_transcript(
    se_tcc = se_igrov_pseudobulk_tcc, em.maxiter = em.maxiter, em.conv = em.conv,
    use_length_normalization = use_length_normalization, ncpu = ncpu
)
saveRDS(se_igrov_pseudobulk_transcript, file.path(result_dir, "se_igrov_pseudobulk_transcript.rds"))

# Preparing the IGROV-1 TCC SE object (k cells subsets)
se_igrov_tcc <- readRDS(file.path(result_dir, "se_igrov_tcc.rds"))
count_matrix <- lapply(k_values, function(k_value) {
    assay(se_igrov_tcc, "counts")[, 1:k_value, drop = FALSE]
})
count_matrix <- do.call(cbind, count_matrix)
se_igrov_k_tcc <- SummarizedExperiment(
    assays = list(counts = count_matrix),
    rowData = rowData(se_igrov_tcc)
)
metadata(se_igrov_k_tcc) <- metadata(se_igrov_tcc)
saveRDS(se_igrov_k_tcc, file.path(result_dir, "se_igrov_k_tcc.rds"))

# Preparing the IGROV-1 pseudobulk TCC SE object (k cells subsets)
se_igrov_k_tcc <- readRDS(file.path(result_dir, "se_igrov_k_tcc.rds"))
cell_labels <- rep(k_values, k_values)
se_igrov_pseudobulk_k_tcc <- pseudobulk_tcc(se_igrov_k_tcc, cell_labels)
saveRDS(se_igrov_pseudobulk_k_tcc, file.path(result_dir, "se_igrov_pseudobulk_k_tcc.rds"))

# Preparing the IGROV-1 pseudobulk gene SE object (k cells subsets)
se_igrov_pseudobulk_k_tcc <- readRDS(file.path(result_dir, "se_igrov_pseudobulk_k_tcc.rds"))
se_igrov_pseudobulk_k_gene <- tcc_to_gene(se_tcc = se_igrov_pseudobulk_k_tcc)
saveRDS(se_igrov_pseudobulk_k_gene, file.path(result_dir, "se_igrov_pseudobulk_k_gene.rds"))

# Preparing the IGROV-1 pseudobulk transcript SE object (k cells subsets)
se_igrov_pseudobulk_k_tcc <- readRDS(file.path(result_dir, "se_igrov_pseudobulk_k_tcc.rds"))
em.maxiter <- 250
em.conv <- 0.01
use_length_normalization <- FALSE
se_igrov_pseudobulk_k_transcript <- tcc_to_transcript(
    se_tcc = se_igrov_pseudobulk_k_tcc, em.maxiter = em.maxiter, em.conv = em.conv,
    use_length_normalization = use_length_normalization, ncpu = ncpu
)
saveRDS(se_igrov_pseudobulk_k_transcript, file.path(result_dir, "se_igrov_pseudobulk_k_transcript.rds"))

# Print session information
sessionInfo()
