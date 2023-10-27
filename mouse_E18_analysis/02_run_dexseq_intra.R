#!/usr/bin/env Rscript

suppressMessages({
    library(tidyverse)
    library(glue)
    library(Isosceles)
    library(DEXSeq)
    library(fitdistrplus)
    library(poolr)
    library(RColorBrewer)
    library(dittoSeq)
    library(Nebulosa)
    require(patchwork)
    library(pdftools)
})

# Set the number of CPUs/threads for the analysis
ncpu <- 1

# Global parameters
window_size <- 30
window_step <- 15
BPPARAM <- MulticoreParam(ncpu)

result_dir <- "02_run_dexseq_intra"
dir.create(result_dir, recursive = TRUE)

# Read the analysis results from previous steps
se_tcc <- readRDS("00_run_isosceles/se_tcc.rds")
se_gene <- readRDS("00_run_isosceles/se_gene.rds")
sce <- readRDS("01_scrnaseq_analysis/sce.rds")
sce_psi <- readRDS("01_scrnaseq_analysis/sce_psi.rds")
pseudotime_matrix <- readRDS("01_scrnaseq_analysis/pseudotime_matrix.rds")
psi_events_list <- readRDS("01_scrnaseq_analysis/psi_events_list.rds")

# Prepare genes of interest for each trajectory
gene_ids_list <- lapply(psi_events_list, function(psi_events) {
    unique(sapply(strsplit(psi_events, ":"), "[", 1))
})

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
    gene_ids <- gene_ids_list[[traj_name]]
    se_window_transcript <- se_window_transcript_list[[traj_name]]
    se_window_psi <- transcript_to_psi(
        se = se_window_transcript[rowData(se_window_transcript)$gene_id %in% gene_ids,],
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
    psi_events <- psi_events_list[[traj_name]]
    se_window_psi <- se_window_psi_list[[traj_name]]
    psi_values <- assay(se_window_psi[psi_events,], "psi")
    return(psi_values)
})
names(psi_values_list) <- colnames(pseudotime_matrix)

gene_counts_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    psi_events <- psi_events_list[[traj_name]]
    se_window_psi <- se_window_psi_list[[traj_name]]
    se_window_gene <- se_window_gene_list[[traj_name]]
    gene_counts <- assay(se_window_gene, "counts")[
        rowData(se_window_psi[psi_events,])$gene_id,
    ]
    return(gene_counts)
})
names(gene_counts_list) <- colnames(pseudotime_matrix)

psi_counts_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    psi_values <- psi_values_list[[traj_name]]
    gene_counts <- gene_counts_list[[traj_name]]
    psi_counts <- psi_values * gene_counts
    # Removing PSI events with redundant count profiles
    gene_ids <- sapply(strsplit(rownames(psi_counts), ":"), "[", 1)
    hash_ids <- apply(psi_counts, 1, rlang::hash)
    row_ids <- paste0(gene_ids, ".", hash_ids)
    row_selector <- !duplicated(row_ids)
    psi_counts <- psi_counts[row_selector,]
    return(psi_counts)
})
names(psi_counts_list) <- colnames(pseudotime_matrix)
saveRDS(psi_counts_list, file.path(result_dir, "psi_counts_list.rds"))

psi_other_counts_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    psi_values <- 1 - psi_values_list[[traj_name]]
    gene_counts <- gene_counts_list[[traj_name]]
    psi_counts <- psi_counts_list[[traj_name]]
    psi_other_counts <- psi_values * gene_counts
    psi_other_counts <- psi_other_counts[rownames(psi_counts),]
    return(psi_other_counts)
})
names(psi_other_counts_list) <- colnames(pseudotime_matrix)
saveRDS(psi_other_counts_list, file.path(result_dir, "psi_other_counts_list.rds"))

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
# (not enough data for the CR cells analysis, so they are excluded)
dexseq_results_list <- lapply(colnames(pseudotime_matrix)[1:5], function(traj_name) {
    psi_counts <- round(as.matrix(psi_counts_list[[traj_name]]))
    psi_other_counts <- round(as.matrix(psi_other_counts_list[[traj_name]]))
    rownames(psi_counts) <- gsub(":", "|", rownames(psi_counts))
    rownames(psi_other_counts) <- gsub(":", "|", rownames(psi_other_counts))
    pseudotime <- window_pseudotime_list[[traj_name]]
    col_data <- data.frame(
        window_name = factor(colnames(psi_counts),
                             levels = colnames(psi_counts)),
        pseudotime = scale(pseudotime)
    )
    dxd <- DEXSeqDataSet(
        countData = psi_counts,
        alternativeCountData = psi_other_counts,
        sampleData = col_data,
        design = ~sample + exon + pseudotime:exon,
        featureID = rownames(psi_counts),
        groupID = rownames(psi_counts)
    )
    dxd <- estimateSizeFactors(dxd)
    dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)
    dxd <- testForDEU(dxd, BPPARAM = BPPARAM)
    dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "pseudotime")
    result_df <- DEXSeqResults(dxd) %>%
        as.data.frame() %>%
        as_tibble() %>%
        dplyr::select(featureID, pvalue, padj) %>%
        transmute(
            psi_event = gsub("\\|", ":", featureID),
            psi_label = prepare_psi_labels(psi_event),
            gene_id = sapply(strsplit(gsub("\\|", ":", featureID), ":"), "[", 1),
            gene_name = rowData(se_gene[gene_id,])$gene_name,
            trajectory = traj_name,
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
    result_df$max_abs_logFC <- max_abs_logFC
    return(result_df)
})
dexseq_results_df <- do.call(rbind, dexseq_results_list)
saveRDS(dexseq_results_df, file.path(result_dir, "dexseq_results_df.rds"))

write.csv(
    filter(dexseq_results_df, fdr <= 0.05, max_abs_logFC >= 1),
    "results/dexseq_results_intra.csv",
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
        perm_psi_counts_list_gene <- lapply(colnames(pseudotime_matrix), function(traj_name) {

            pseudotime <- pseudotime_matrix[, traj_name]
            pseudotime <- pseudotime[!is.na(pseudotime)]
            se_tcc_traj <- se_tcc[, names(pseudotime)]
            se_gene_traj <- se_gene[, names(pseudotime)]
            gene_ids <- gene_ids_list[[traj_name]]
            se_window_gene <- se_window_gene_list[[traj_name]]

            tcc_row_selector <- rowData(se_tcc_traj)$gene_id == gene_id
            tcc_tx_selector <- metadata(se_tcc_traj)$transcript_df$gene_id == gene_id
            perm_se_tcc <- se_tcc_traj[tcc_row_selector,]
            metadata(perm_se_tcc)$compatibility_matrix <-
                metadata(perm_se_tcc)$compatibility_matrix[tcc_row_selector, tcc_tx_selector]
            metadata(perm_se_tcc)$transcript_df <-
                metadata(perm_se_tcc)$transcript_df[tcc_tx_selector,]
            metadata(perm_se_tcc)$transcript_exon_granges_list <-
                metadata(perm_se_tcc)$transcript_exon_granges_list[tcc_tx_selector,]
            non_zero_indices <- which(assay(se_gene_traj[gene_id,], "counts")[1,] != 0)
            if (length(non_zero_indices) > 1) {
                cell_indices <- seq(ncol(perm_se_tcc))
                cell_indices[non_zero_indices] <- sample(cell_indices[non_zero_indices])
                perm_se_tcc <- perm_se_tcc[, cell_indices]
            }

            perm_se_window_tcc <- pseudotime_tcc(
                se_tcc = perm_se_tcc,
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

# Calculate p-values for PSI events of the Celf2 gene from permuted PSI counts

filtered_results_df <- filter(dexseq_results_df, gene_name == "Celf2")
perm_pvalues <- sapply(seq(nrow(filtered_results_df)), function(i) {
    psi_event <- filtered_results_df$psi_event[i]
    gene_id <- filtered_results_df$gene_id[i]
    traj_name <- filtered_results_df$trajectory[i]
    psi_counts <- psi_counts_list[[traj_name]]
    window_number <- ncol(psi_counts)
    eff_test_number <- 1 + (window_number - 1) * window_step / window_size
    psi_event_counts <- psi_counts[psi_event,]
    perm_psi_event_counts <- sapply(seq(n_perm), function(i) {
        perm_psi_counts_list[[gene_id]][[i]][[traj_name]][psi_event,]
    })
    psi_event_pvalues <- sapply(names(psi_event_counts), function(window_id) {
        if (all(round(perm_psi_event_counts[window_id,]) == 0)) return(1)
        fit_nb <- fitdist(round(perm_psi_event_counts[window_id,]), "nbinom")
        hi_val <- pnbinom(psi_event_counts[window_id],
                          mu = fit_nb$estimate[2], size = fit_nb$estimate[1],
                          lower.tail = FALSE)
        lo_val <- pnbinom(psi_event_counts[window_id],
                          mu = fit_nb$estimate[2], size = fit_nb$estimate[1],
                          lower.tail = TRUE)
        return(c(lo_val, hi_val))
    })
    rownames(psi_event_pvalues) <- c("low", "high")
    psi_event_pvalue <- 2 * min(c(
        fisher(c(psi_event_pvalues["low", 1],
                 psi_event_pvalues["high", ncol(psi_event_pvalues)]))$p,
        fisher(c(psi_event_pvalues["high", 1],
                 psi_event_pvalues["low", ncol(psi_event_pvalues)]))$p
    ))
    return(psi_event_pvalue)
})
filtered_results_df$perm_pvalue <- perm_pvalues

# Celf2:0c9e:A5 permutation p-value (0.000167)
filtered_results_df %>%
    filter(psi_label == "Celf2:0c9e:A5", trajectory == "glut_1") %>%
    pull(perm_pvalue)

# Celf2:fc81:A3 permutation p-value (0.206)
filtered_results_df %>%
    filter(psi_label == "Celf2:fc81:A3", trajectory == "glut_1") %>%
    pull(perm_pvalue)

# Combined p-value (0.000388)
fisher(c(A5 = 0.000167, A3 = 0.206))$p

################################################################################

# Create gene & PSI expression UMAP plots for significant results

plot_umap_expression <- function(psi_event, psi_label, pal = NULL) {
    gene_id <- strsplit(psi_event, ":")[[1]][1]
    p1 <- plot_density(sce_psi, psi_event,
                       slot = "psi", size = 1.5) +
        labs(title = "PSI values density") +
        theme(legend.position = "right",
              plot.title = element_text(size = 13),
              legend.title = element_text(size = 11))
    p2 <- plot_density(sce, gene_id,
                       size = 1.5) +
        labs(title = "Gene expression density") +
        theme(legend.position = "right",
              plot.title = element_text(size = 13),
              legend.title = element_text(size = 11))
    if (!is.null(pal)) {
        p1 <- p1 +
            scale_color_gradientn(
                colours = brewer.pal(n = 7, name = pal)
            )
        p2 <- p2 +
            scale_color_gradientn(
                colours = brewer.pal(n = 7, name = pal)
            )
    }
    p3 <- dittoDimPlot(sce_psi,
                       psi_event,
                       reduction.use = "UMAP",
                       assay = "psi",
                       size = 1.5, main = "PSI values",
                       legend.title = "PSI",
                       order = "increasing")
    p4 <- dittoDimPlot(sce, gene_id,
                       reduction.use = "UMAP",
                       size = 1.5, main = "Gene expression",
                       legend.title = "Logcounts",
                       order = "increasing")

    patchwork <- (p1 + p2) / (p3 + p4)
    patchwork + plot_annotation(title = psi_label, subtitle = psi_event)
}

# Expression UMAP plots for selected PSI events of the Celf2 gene
umap_plot <- plot_umap_expression(
    "ENSMUSG00000002107:chr2:6560659-6560670:-:A5", "Celf2:0c9e:A5",
    pal = "Reds"
)
ggsave("results/umap_psi_Celf2_A5.pdf", umap_plot, height = 5, width = 8)

umap_plot <- plot_umap_expression(
    "ENSMUSG00000002107:chr2:6553965-6553982:-:A3", "Celf2:fc81:A3",
    pal = "Blues"
)
ggsave("results/umap_psi_Celf2_A3.pdf", umap_plot, height = 5, width = 8)

umap_plot <- plot_umap_expression(
    "ENSMUSG00000002107:chr2:6546780-6547041:-:RI", "Celf2:823a:RI",
    pal = "Greens"
)
ggsave("results/umap_psi_Celf2_RI.pdf", umap_plot, height = 5, width = 8)

# Expression UMAP plots for all identified significant PSI events
filtered_results_df <- filter(dexseq_results_df, fdr <= 0.05, max_abs_logFC >= 1)
psi_label_to_event <- setNames(
    filtered_results_df$psi_event, filtered_results_df$psi_label
)
psi_labels <- sort(unique(filtered_results_df$psi_label))
dir.create("results/temp", recursive = TRUE)
pdf_files <- sapply(seq_along(psi_labels), function(i) {
    psi_label <- psi_labels[i]
    psi_event <- unname(psi_label_to_event[psi_label])
    umap_plot <- plot_umap_expression(psi_event, psi_label)
    pdf_file <- glue("results/temp/{i}.pdf")
    ggsave(pdf_file, umap_plot, height = 5, width = 8)
    return(pdf_file)
})
pdf_combine(pdf_files, output = "results/umap_psi_intra.pdf")
unlink("results/temp", recursive = TRUE)

################################################################################

# Create PSI counts trajectory plots for significant results

trajectory_names <- c(
    glut_1 = "Glutamatergic trajectory 1",
    glut_2 = "Glutamatergic trajectory 2",
    gaba = "GABAergic trajectory",
    rad_glia = "Radial glia",
    cyc_rad_glia = "Cycling radial glia",
    cr = "Cajal Retzius cells"
)

plot_psi_single_traj <- function(psi_event, psi_label, traj_name, plot_title = TRUE) {
    gene_id <- strsplit(psi_event, ":")[[1]][1]
    pseudotime <- window_pseudotime_list[[traj_name]]
    psi_df <- data.frame(
        pseudotime = pseudotime,
        psi_counts = plot_psi_counts_list[[traj_name]][psi_event,]
    )
    psi_perm_matrix <- sapply(seq(n_perm), function(i) {
        perm_matrix <- perm_psi_counts_list[[gene_id]][[i]][[traj_name]]
        if (psi_event %in% rownames(perm_matrix)) {
            return(perm_matrix[psi_event,])
        } else {
            return(setNames(rep(0, ncol(perm_matrix)), colnames(perm_matrix)))
        }
    })
    psi_perm_df <- data.frame(
        pseudotime = rep(pseudotime, ncol(psi_perm_matrix)),
        psi_counts = as.vector(psi_perm_matrix),
        permutation = rep(seq(ncol(psi_perm_matrix)), each = length(pseudotime))
    )
    ggplot() +
        geom_line(
            data = psi_perm_df,
            col = "grey", alpha = 0.5,
            mapping = aes(x = pseudotime, y = psi_counts, group = permutation)
        ) +
        geom_line(
            data = psi_df,
            col = "red", size = 1,
            mapping = aes(x = pseudotime, y = psi_counts)
        ) +
        labs(
            title = ifelse(plot_title, psi_label, ""),
            subtitle = ifelse(plot_title, psi_event, trajectory_names[traj_name]),
            x = "Mean window pseudotime",
            y = "PSI counts"
        ) +
        theme_bw()
}

plot_psi_all_traj <- function(psi_event, psi_label) {
    p1 <- plot_psi_single_traj(psi_event, psi_label, "glut_1", plot_title = FALSE)
    p2 <- plot_psi_single_traj(psi_event, psi_label, "glut_2", plot_title = FALSE)
    p3 <- plot_psi_single_traj(psi_event, psi_label, "gaba", plot_title = FALSE)
    p4 <- plot_psi_single_traj(psi_event, psi_label, "rad_glia", plot_title = FALSE)
    p5 <- plot_psi_single_traj(psi_event, psi_label, "cyc_rad_glia", plot_title = FALSE)
    p6 <- plot_psi_single_traj(psi_event, psi_label, "cr", plot_title = FALSE)
    p1 + p4 + p2 + p5 + p3 + p6 +
        plot_layout(nrow = 3, widths = c(2, 1)) +
        plot_annotation(title = psi_label, subtitle = psi_event)
}

# PSI counts trajectory plots for selected PSI events of the Celf2 gene
## Glutamatergic trajectory 1
p1 <- plot_psi_single_traj("ENSMUSG00000002107:chr2:6560659-6560670:-:A5",
                           "Celf2:0c9e:A5", "glut_1")
p2 <- plot_psi_single_traj("ENSMUSG00000002107:chr2:6553965-6553982:-:A3",
                           "Celf2:fc81:A3", "glut_1")
p3 <- plot_psi_single_traj("ENSMUSG00000002107:chr2:6546780-6547041:-:RI",
                           "Celf2:823a:RI", "glut_1")
traj_plot <- p1 / p2 / p3
ggsave("results/trajectory_plot_intra_Celf2_Glut_1.pdf", traj_plot,
       height = 8, width = 8)
## All trajectories
traj_plot <- plot_psi_all_traj("ENSMUSG00000002107:chr2:6560659-6560670:-:A5",
                               "Celf2:0c9e:A5")
ggsave("results/trajectory_plot_intra_Celf2_A5.pdf", traj_plot,
       height = 8, width = 8)
traj_plot <- plot_psi_all_traj("ENSMUSG00000002107:chr2:6553965-6553982:-:A3",
                               "Celf2:fc81:A3")
ggsave("results/trajectory_plot_intra_Celf2_A3.pdf", traj_plot,
       height = 8, width = 8)
traj_plot <- plot_psi_all_traj("ENSMUSG00000002107:chr2:6546780-6547041:-:RI",
                               "Celf2:823a:RI")
ggsave("results/trajectory_plot_intra_Celf2_RI.pdf", traj_plot,
       height = 8, width = 8)

# PSI counts trajectory plots for all identified significant PSI events
filtered_results_df <- filter(dexseq_results_df, fdr <= 0.05, max_abs_logFC >= 1)
psi_labels <- sort(unique(filtered_results_df$psi_label))
psi_label_to_event <- setNames(
    filtered_results_df$psi_event, filtered_results_df$psi_label
)
dir.create("results/temp", recursive = TRUE)
pdf_files <- sapply(seq_along(psi_labels), function(i) {
    psi_label <- psi_labels[i]
    psi_event <- unname(psi_label_to_event[psi_label])
    traj_plot <- plot_psi_all_traj(psi_event, psi_label)
    pdf_file <- glue("results/temp/{i}.pdf")
    ggsave(pdf_file, traj_plot, height = 8, width = 8)
    return(pdf_file)
})
pdf_combine(pdf_files, output = "results/trajectory_plot_intra.pdf")
unlink("results/temp", recursive = TRUE)

################################################################################

# Print session information
sessionInfo()
