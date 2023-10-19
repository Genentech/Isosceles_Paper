#!/usr/bin/env Rscript

suppressMessages({
    library(tidyverse)
    library(glue)
    library(Isosceles)
    library(DEXSeq)
})

# Set the number of CPUs/threads for the analysis
ncpu <- 1

# Global parameters
window_size <- 30
window_step <- 15
BPPARAM <- MulticoreParam(ncpu)

result_dir <- "03_run_dexseq_inter"
dir.create(result_dir, recursive = TRUE)

# Read the analysis results from previous steps
se_tcc <- readRDS("00_run_isosceles/se_tcc.rds")
se_gene <- readRDS("00_run_isosceles/se_gene.rds")
pseudotime_matrix <- readRDS("01_scrnaseq_analysis/pseudotime_matrix.rds")
psi_events_list <- readRDS("01_scrnaseq_analysis/psi_events_list.rds")

# Prepare shared PSI regions and genes of interest for each trajectory
shared_psi_events <- table(unlist(psi_events_list))
shared_psi_events <- names(shared_psi_events[shared_psi_events > 1])
shared_gene_ids <- unique(sapply(strsplit(shared_psi_events, ":"), "[", 1))

################################################################################

# Pseudotime window analysis

se_window_tcc_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    pseudotime <- pseudotime_matrix[, traj_name]
    pseudotime <- pseudotime[!is.na(pseudotime)]
    se_tcc_traj <- se_tcc[, names(pseudotime)]
    se_window_tcc <- pseudotime_tcc(
        se_tcc = se_tcc_traj,
        pseudotime = pseudotime,
        trim = 0,
        window_size = window_size,
        window_step = window_step
    )
    return(se_window_tcc)
})
names(se_window_tcc_list) <- colnames(pseudotime_matrix)

se_window_gene_list <- lapply(se_window_tcc_list, function(se_window_tcc) {
    se_window_gene <- tcc_to_gene(
        se_tcc = se_window_tcc
    )
    return(se_window_gene)
})

se_window_transcript_list <- lapply(se_window_tcc_list, function(se_window_tcc) {
    se_window_transcript <- tcc_to_transcript(
        se_tcc = se_window_tcc,
        em.maxiter = 250, em.conv = 0.01,
        use_length_normalization = FALSE, ncpu = ncpu
    )
    return(se_window_transcript)
})

se_window_psi_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    se_window_transcript <- se_window_transcript_list[[traj_name]]
    se_window_psi <- transcript_to_psi(
        se = se_window_transcript[rowData(se_window_transcript)$gene_id %in% shared_gene_ids,],
        ncpu = ncpu
    )
    return(se_window_psi)
})
names(se_window_psi_list) <- colnames(pseudotime_matrix)

window_pseudotime_list <- lapply(se_window_tcc_list, function(se_window_tcc) {
    return(se_window_tcc$pseudotime)
})
saveRDS(window_pseudotime_list, file.path(result_dir, "window_pseudotime_list.rds"))

################################################################################

# Prepare PSI count data

psi_values_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    se_window_psi <- se_window_psi_list[[traj_name]]
    psi_values <- as.matrix(assay(se_window_psi, "psi"))
    psi_values <- sapply(colnames(psi_values), function(window_id) {
        psi_vector <- psi_values[, window_id]
        psi_vector <- psi_vector[shared_psi_events]
        psi_vector[is.na(psi_vector)] <- 0
        names(psi_vector) <- shared_psi_events
        return(psi_vector)
    })
    return(psi_values)
})
names(psi_values_list) <- colnames(pseudotime_matrix)

gene_counts_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    se_window_gene <- se_window_gene_list[[traj_name]]
    gene_counts <- as.matrix(assay(se_window_gene, "counts")[
        sapply(strsplit(shared_psi_events, ":"), "[", 1),
    ])
    return(gene_counts)
})
names(gene_counts_list) <- colnames(pseudotime_matrix)

psi_counts_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    psi_values <- psi_values_list[[traj_name]]
    gene_counts <- gene_counts_list[[traj_name]]
    psi_counts <- psi_values * gene_counts
    colnames(psi_counts) <- paste0(traj_name, ".", colnames(psi_counts))
    return(psi_counts)
})
names(psi_counts_list) <- colnames(pseudotime_matrix)
saveRDS(psi_counts_list, file.path(result_dir, "psi_counts_list.rds"))

psi_other_counts_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    psi_values <- 1 - psi_values_list[[traj_name]]
    gene_counts <- gene_counts_list[[traj_name]]
    psi_other_counts <- psi_values * gene_counts
    colnames(psi_other_counts) <- paste0(traj_name, ".", colnames(psi_other_counts))
    return(psi_other_counts)
})
names(psi_other_counts_list) <- colnames(pseudotime_matrix)

psi_counts <- do.call(cbind, psi_counts_list)
psi_other_counts <- do.call(cbind, psi_other_counts_list)
# Remove PSI events with redundant count profiles
psi_gene_ids <- sapply(strsplit(rownames(psi_counts), ":"), "[", 1)
psi_hash_ids <- apply(psi_counts, 1, rlang::hash)
psi_row_ids <- paste0(psi_gene_ids, ".", psi_hash_ids)
row_selector <- !duplicated(psi_row_ids)
psi_counts <- psi_counts[row_selector,]
psi_other_counts <- psi_other_counts[row_selector,]
saveRDS(psi_counts, file.path(result_dir, "psi_counts.rds"))
saveRDS(psi_other_counts, file.path(result_dir, "psi_other_counts.rds"))

################################################################################

# Helper function for the DEXSeq analysis
prepare_psi_labels <- function(psi_events) {
    gene_id_to_name <- setNames(rowData(se_gene)$gene_name,
                                rowData(se_gene)$gene_id)
    gene_ids <- sapply(strsplit(psi_events, ":"), "[", 1)
    gene_names <- unname(gene_id_to_name[gene_ids])
    region_types <- sapply(strsplit(psi_events, ":"), "[", 5)
    hash_ids <- unname(sapply(psi_events, rlang::hash))
    psi_labels <- paste0(gene_names, ":",
                         substring(hash_ids, 1, 4), ":",
                         region_types)
    return(psi_labels)
}

################################################################################

# Run DEXSeq
dexseq_psi_counts <- round(psi_counts)
dexseq_psi_other_counts <- round(psi_other_counts)
rownames(dexseq_psi_counts) <- gsub(":", "|", rownames(dexseq_psi_counts))
rownames(dexseq_psi_other_counts) <- gsub(":", "|", rownames(dexseq_psi_other_counts))
dexseq_pseudotime <- do.call(c, lapply(colnames(pseudotime_matrix), function(traj_name) {
    pseudotime <- window_pseudotime_list[[traj_name]]
    pseudotime <- scale(pseudotime)
    return(pseudotime)
}))
col_data <- data.frame(
    window_name = factor(colnames(dexseq_psi_counts),
                         levels = colnames(dexseq_psi_counts)),
    trajectory = factor(sapply(strsplit(colnames(dexseq_psi_counts), "\\."), "[", 1),
                        levels = colnames(pseudotime_matrix)),
    pseudotime = dexseq_pseudotime
)
dxd <- DEXSeqDataSet(
    countData = dexseq_psi_counts,
    alternativeCountData = dexseq_psi_other_counts,
    sampleData = col_data,
    design = ~sample + exon + pseudotime:exon + trajectory:exon,
    featureID = rownames(dexseq_psi_counts),
    groupID = rownames(dexseq_psi_counts)
)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)
dxd <- testForDEU(dxd, BPPARAM = BPPARAM)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "pseudotime")
dexseq_results_df <- DEXSeqResults(dxd) %>%
    as.data.frame() %>%
    as_tibble() %>%
    dplyr::select(featureID, pvalue, padj) %>%
    transmute(
        psi_event = gsub("\\|", ":", featureID),
        psi_label = prepare_psi_labels(psi_event),
        gene_id = sapply(strsplit(gsub("\\|", ":", featureID), ":"), "[", 1),
        gene_name = rowData(se_gene[gene_id,])$gene_name,
        pvalue = pvalue,
        fdr = padj
    )
logFC_values <- DEXSeqResults(dxd) %>%
    as.data.frame()
logFC_values <- logFC_values[, grepl("^log2fold_", colnames(logFC_values))]
logFC_values <- as.matrix(logFC_values)
max_abs_logFC <- apply(logFC_values, 1, function(x) {max(abs(x))})
max_abs_logFC[is.na(max_abs_logFC)] <- 0
names(max_abs_logFC) <- NULL
dexseq_results_df$max_abs_logFC <- max_abs_logFC
saveRDS(dexseq_results_df, file.path(result_dir, "dexseq_results_df.rds"))

write.csv(
    filter(dexseq_results_df, fdr <= 0.05, max_abs_logFC >= 1),
    "results/dexseq_results_inter.csv",
    row.names = FALSE
)

################################################################################

# Recalculate PSI count data for trajectory plots

filtered_results_df <- filter(dexseq_results_df, fdr <= 0.05, max_abs_logFC >= 1)
filtered_gene_ids <- unique(filtered_results_df$gene_id)
plot_se_window_psi_list <- lapply(se_window_transcript_list, function(se) {
    se_window_psi <- transcript_to_psi(
        se = se[rowData(se)$gene_id %in% filtered_gene_ids,],
        ncpu = ncpu
    )
    return(se_window_psi)
})
filtered_psi_events <- unique(unlist(lapply(plot_se_window_psi_list, rownames)))

plot_psi_counts_list <- lapply(names(plot_se_window_psi_list), function(traj_name) {
    se_window_psi <- plot_se_window_psi_list[[traj_name]]
    se_window_gene <- se_window_gene_list[[traj_name]]
    psi_values <- as.matrix(assay(se_window_psi, "psi"))
    psi_values <- sapply(colnames(psi_values), function(window_id) {
        psi_vector <- psi_values[, window_id]
        psi_vector <- psi_vector[filtered_psi_events]
        psi_vector[is.na(psi_vector)] <- 0
        names(psi_vector) <- filtered_psi_events
        return(psi_vector)
    })
    gene_counts <- as.matrix(assay(se_window_gene, "counts")[
        sapply(strsplit(filtered_psi_events, ":"), "[", 1),
    ])
    psi_counts <- psi_values * gene_counts
    return(psi_counts)
})
names(plot_psi_counts_list) <- names(plot_se_window_psi_list)
saveRDS(plot_psi_counts_list, file.path(result_dir, "plot_psi_counts_list.rds"))

################################################################################

# Prepare permuted PSI count data

filtered_results_df <- filter(dexseq_results_df, fdr <= 0.05)
filtered_results_df <- filtered_results_df[!duplicated(filtered_results_df$psi_event),]
filtered_psi_events <- filtered_results_df$psi_event
filtered_gene_ids <- unique(filtered_results_df$gene_id)
n_perm <- 100

set.seed(42)
perm_psi_counts_list <- lapply(filtered_gene_ids, function(gene_id) {
    lapply(seq(n_perm), function(i) {
        orig_se_gene <- se_gene[, rownames(pseudotime_matrix)]
        perm_se_tcc <- se_tcc[, rownames(pseudotime_matrix)]
        tcc_row_selector <- rowData(perm_se_tcc)$gene_id == gene_id
        tcc_tx_selector <- metadata(perm_se_tcc)$transcript_df$gene_id == gene_id
        perm_se_tcc <- perm_se_tcc[tcc_row_selector,]
        metadata(perm_se_tcc)$compatibility_matrix <-
            metadata(perm_se_tcc)$compatibility_matrix[tcc_row_selector, tcc_tx_selector]
        metadata(perm_se_tcc)$transcript_df <-
            metadata(perm_se_tcc)$transcript_df[tcc_tx_selector,]
        metadata(perm_se_tcc)$transcript_exon_granges_list <-
            metadata(perm_se_tcc)$transcript_exon_granges_list[tcc_tx_selector,]
        non_zero_indices <- which(assay(orig_se_gene[gene_id,], "counts")[1,] != 0)
        cell_indices <- seq(ncol(perm_se_tcc))
        cell_indices[non_zero_indices] <- sample(cell_indices[non_zero_indices])
        cell_ids <- colnames(perm_se_tcc)
        perm_se_tcc <- perm_se_tcc[, cell_indices]
        colnames(perm_se_tcc) <- cell_ids

        perm_psi_counts_list_gene <- lapply(colnames(pseudotime_matrix), function(traj_name) {

            pseudotime <- pseudotime_matrix[, traj_name]
            pseudotime <- pseudotime[!is.na(pseudotime)]
            perm_se_tcc_traj <- perm_se_tcc[, names(pseudotime)]
            se_window_gene <- se_window_gene_list[[traj_name]]

            perm_se_window_tcc <- pseudotime_tcc(
                se_tcc = perm_se_tcc_traj,
                pseudotime = pseudotime,
                trim = 0,
                window_size = window_size,
                window_step = window_step
            )
            perm_se_window_transcript <- tcc_to_transcript(
                se_tcc = perm_se_window_tcc,
                em.maxiter = 250, em.conv = 0.01,
                use_length_normalization = FALSE, ncpu = ncpu
            )
            perm_se_window_psi <- transcript_to_psi(
                perm_se_window_transcript,
                ncpu = ncpu
            )

            perm_psi_values <- assay(perm_se_window_psi, "psi")
            orig_gene_counts <- assay(se_window_gene, "counts")[rowData(perm_se_window_psi)$gene_id,]

            return(perm_psi_values * orig_gene_counts)
        })
        names(perm_psi_counts_list_gene) <- colnames(pseudotime_matrix)
        return(perm_psi_counts_list_gene)
    })
})
names(perm_psi_counts_list) <- filtered_gene_ids
saveRDS(perm_psi_counts_list, file.path(result_dir, "perm_psi_counts_list.rds"))

################################################################################

# Print session information
sessionInfo()
