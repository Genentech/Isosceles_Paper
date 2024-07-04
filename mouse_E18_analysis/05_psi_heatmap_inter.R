#!/usr/bin/env Rscript

suppressMessages({
    library(tidyverse)
    library(glue)
    library(Isosceles)
    library(RColorBrewer)
    library(pheatmap)
})

# Set the number of CPUs/threads for the analysis
ncpu <- 1

# Global parameters
heatmap_window_sizes <- c(
    glut_1 = 100, glut_2 = 100, gaba = 100,
    rad_glia = 50, cyc_rad_glia = 50, cr = 50
)
heatmap_window_steps <- c(
    glut_1 = 3, glut_2 = 3, gaba = 3,
    rad_glia = 3, cyc_rad_glia = 3, cr = 3
)
heatmap_n_perm <- 100

result_dir <- "05_psi_heatmap_inter"
dir.create(result_dir, recursive = TRUE)

# Read the analysis results from previous steps
se_tcc <- readRDS("00_run_isosceles/se_tcc.rds")
se_gene <- readRDS("00_run_isosceles/se_gene.rds")
pseudotime_matrix <- readRDS("01_scrnaseq_analysis/pseudotime_matrix.rds")
dexseq_results_df <- readRDS("03_run_dexseq_inter/dexseq_results_df.rds")

################################################################################

# Filter DEXSeq results
filtered_results_df <- filter(dexseq_results_df, fdr <= 0.05)
filtered_results_df <- filtered_results_df[!duplicated(filtered_results_df$psi_event),]
heatmap_psi_events <- filtered_results_df$psi_event
heatmap_gene_ids <- unique(filtered_results_df$gene_id)

################################################################################

# Pseudotime window analysis

tcc_row_selector <- rowData(se_tcc)$gene_id %in% heatmap_gene_ids
tcc_tx_selector <- metadata(se_tcc)$transcript_df$gene_id %in% heatmap_gene_ids
heatmap_se_tcc <- se_tcc[tcc_row_selector,]
metadata(heatmap_se_tcc)$compatibility_matrix <-
    metadata(heatmap_se_tcc)$compatibility_matrix[tcc_row_selector, tcc_tx_selector]
metadata(heatmap_se_tcc)$transcript_df <-
    metadata(heatmap_se_tcc)$transcript_df[tcc_tx_selector,]
metadata(heatmap_se_tcc)$transcript_exon_granges_list <-
    metadata(heatmap_se_tcc)$transcript_exon_granges_list[tcc_tx_selector,]

heatmap_se_window_tcc_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    pseudotime <- pseudotime_matrix[, traj_name]
    pseudotime <- pseudotime[!is.na(pseudotime)]
    se_tcc_traj <- heatmap_se_tcc[, names(pseudotime)]
    se_window_tcc <- pseudotime_tcc(
        se_tcc = se_tcc_traj,
        pseudotime = pseudotime,
        trim = 0,
        window_size = heatmap_window_sizes[traj_name],
        window_step = heatmap_window_steps[traj_name]
    )
    return(se_window_tcc)
})
names(heatmap_se_window_tcc_list) <- colnames(pseudotime_matrix)

heatmap_se_window_gene_list <- lapply(heatmap_se_window_tcc_list, function(se_window_tcc) {
    se_window_gene <- tcc_to_gene(
        se_tcc = se_window_tcc
    )
    return(se_window_gene)
})

heatmap_se_window_transcript_list <- lapply(heatmap_se_window_tcc_list, function(se_window_tcc) {
    se_window_transcript <- tcc_to_transcript(
        se_tcc = se_window_tcc,
        em.maxiter = 250, em.conv = 0.01,
        use_length_normalization = FALSE, ncpu = ncpu
    )
    return(se_window_transcript)
})

heatmap_se_window_psi_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    se_window_transcript <- heatmap_se_window_transcript_list[[traj_name]]
    se_window_psi <- transcript_to_psi(
        se = se_window_transcript,
        ncpu = ncpu
    )
    return(se_window_psi)
})
names(heatmap_se_window_psi_list) <- colnames(pseudotime_matrix)

heatmap_window_pseudotime_list <- lapply(heatmap_se_window_tcc_list, function(se_window_tcc) {
    return(se_window_tcc$pseudotime)
})
saveRDS(heatmap_window_pseudotime_list, file.path(result_dir, "heatmap_window_pseudotime_list.rds"))

################################################################################

# Prepare PSI count data

heatmap_psi_values_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    se_window_psi <- heatmap_se_window_psi_list[[traj_name]]
    psi_values <- as.matrix(assay(se_window_psi, "psi"))
    psi_values <- sapply(colnames(psi_values), function(window_id) {
        psi_vector <- psi_values[, window_id]
        psi_vector <- psi_vector[heatmap_psi_events]
        psi_vector[is.na(psi_vector)] <- 0
        names(psi_vector) <- heatmap_psi_events
        return(psi_vector)
    })
    return(psi_values)
})
names(heatmap_psi_values_list) <- colnames(pseudotime_matrix)

heatmap_gene_counts_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    se_window_gene <- heatmap_se_window_gene_list[[traj_name]]
    gene_counts <- as.matrix(assay(se_window_gene, "counts")[
        sapply(strsplit(heatmap_psi_events, ":"), "[", 1),
    ])
    return(gene_counts)
})
names(heatmap_gene_counts_list) <- colnames(pseudotime_matrix)

heatmap_psi_counts_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    psi_values <- heatmap_psi_values_list[[traj_name]]
    gene_counts <- heatmap_gene_counts_list[[traj_name]]
    psi_counts <- psi_values * gene_counts
    return(psi_counts)
})
names(heatmap_psi_counts_list) <- colnames(pseudotime_matrix)
saveRDS(heatmap_psi_counts_list, file.path(result_dir, "heatmap_psi_counts_list.rds"))

################################################################################

# Prepare permuted PSI count data

set.seed(42)
heatmap_perm_psi_counts_list <- lapply(heatmap_gene_ids, function(gene_id) {
    lapply(seq(heatmap_n_perm), function(i) {
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
            se_window_gene <- heatmap_se_window_gene_list[[traj_name]]
            perm_se_window_tcc <- pseudotime_tcc(
                se_tcc = perm_se_tcc_traj,
                pseudotime = pseudotime,
                trim = 0,
                window_size = heatmap_window_sizes[traj_name],
                window_step = heatmap_window_steps[traj_name]
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
names(heatmap_perm_psi_counts_list) <- heatmap_gene_ids
saveRDS(heatmap_perm_psi_counts_list, file.path(result_dir, "heatmap_perm_psi_counts_list.rds"))

################################################################################

# Prepare heatmap data

heatmap_psi_counts <- do.call(cbind, lapply(names(heatmap_psi_counts_list), function(traj_name) {
    psi_counts <- heatmap_psi_counts_list[[traj_name]]
    colnames(psi_counts) <- paste0(traj_name, ".", colnames(psi_counts))
    return(as.matrix(psi_counts))
}))

heatmap_avg_perm_psi_counts_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    avg_perm_matrix <- t(sapply(heatmap_psi_events, function(psi_event) {
        gene_id <- strsplit(psi_event, ":")[[1]][1]
        perm_psi_event_counts <- sapply(seq(heatmap_n_perm), function(i) {
            perm_matrix <- heatmap_perm_psi_counts_list[[gene_id]][[i]][[traj_name]]
            if (psi_event %in% rownames(perm_matrix)) {
                return(perm_matrix[psi_event,])
            } else {
                return(setNames(rep(0, ncol(perm_matrix)), colnames(perm_matrix)))
            }
        })
        avg_psi_event_counts <- apply(perm_psi_event_counts, 1, mean)
        names(avg_psi_event_counts) <- paste0(
            traj_name, ".", rownames(perm_psi_event_counts)
        )
        return(avg_psi_event_counts)
    }))
})
names(heatmap_avg_perm_psi_counts_list) <- colnames(pseudotime_matrix)
heatmap_avg_perm_psi_counts <- do.call(cbind, heatmap_avg_perm_psi_counts_list)

heatmap_values <- heatmap_psi_counts / heatmap_avg_perm_psi_counts
heatmap_values <- log2(heatmap_values + 0.1)
heatmap_values[is.nan(heatmap_values)] <- 0
saveRDS(heatmap_values, file.path(result_dir, "heatmap_values.rds"))

################################################################################

# Create the heatmap

traj_colors <- list(
    glut_1 = brewer.pal(n = 7, name = "YlGnBu")[2:6],
    glut_2 = brewer.pal(n = 7, name = "YlGn")[2:6],
    gaba = brewer.pal(n = 7, name = "Reds")[2:6],
    rad_glia = c("white", "deeppink"),
    cyc_rad_glia = c("white", "purple"),
    cr = c("white", "gold")
)

is_correlated_row <- function(mat) {
    cor_mat <- cor(t(mat))
    cor_mat[!lower.tri(cor_mat)] <- 0
    apply(cor_mat, 2, function(x) any(abs(x) > 0.99, na.rm = TRUE))
}

filtered_results_df <- filter(dexseq_results_df, fdr <= 0.05, max_abs_logFC >= 1)
heatmap_values <- heatmap_values[
    rownames(heatmap_values) %in% filtered_results_df$psi_event,
]
heatmap_values <- heatmap_values[!is_correlated_row(heatmap_values),]
psi_event_to_label <- setNames(
    filtered_results_df$psi_label, filtered_results_df$psi_event
)
rownames(heatmap_values) <- unname(psi_event_to_label[rownames(heatmap_values)])

col_ann_df <- lapply(names(heatmap_se_window_tcc_list), function(traj_name) {
    se <- heatmap_se_window_tcc_list[[traj_name]]
    pseudotime <- se$pseudotime
    names(pseudotime) <- paste0(traj_name, ".", colnames(se))
    pseudotime <- pseudotime[colnames(heatmap_values)]
    names(pseudotime) <- colnames(heatmap_values)
    return(pseudotime)

})
names(col_ann_df) <- names(heatmap_se_window_tcc_list)
col_ann_df <- as.data.frame(do.call(cbind, col_ann_df))

psi_heatmap <- pheatmap(
    heatmap_values,
    cutree_rows = 2,
    breaks = seq(-3, 3, length.out = 101),
    cluster_rows = TRUE, cluster_cols = FALSE,
    show_rownames = TRUE, show_colnames = FALSE,
    annotation_col = col_ann_df,
    annotation_colors = traj_colors,
    gaps_col = cumsum(sapply(heatmap_window_pseudotime_list, length)),
    legend = TRUE, annotation_legend = FALSE,
    annotation_names_col = FALSE,
    fontsize_row = 5, treeheight_row = 0,
    scale = "row"
)
ggsave("results/psi_heatmap_inter.pdf",
       psi_heatmap, height = 3.1, width = 7.1)

################################################################################

# Print session information
sessionInfo()
