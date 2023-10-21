#!/usr/bin/env Rscript

suppressMessages({
    library(tidyverse)
    library(glue)
    library(scran)
    library(BiocParallel)
    library(dittoSeq)
    library(Nebulosa)
    require(patchwork)
    library(pdftools)
})

# Set the number of CPUs/threads for the analysis
ncpu <- 1

# Global parameters
BPPARAM <- MulticoreParam(ncpu)

result_dir <- "07_isoswitch_analysis"
dir.create(result_dir, recursive = TRUE)

# Read the analysis results from previous steps
sce <- readRDS("01_scrnaseq_analysis/sce.rds")
sce_transcript <- readRDS("01_scrnaseq_analysis/sce_transcript.rds")

################################################################################

# Identify marker transcripts between each pair of clusters

sce_transcript_marker <- sce_transcript
sce_transcript_marker <- sce_transcript_marker[
    !is.na(rowData(sce_transcript_marker)$compatible_tx),
]
sce_transcript_marker <- sce_transcript_marker[
    rowMeans(assay(sce_transcript_marker, "tpm")) >= 10,
]
sce_transcript_marker$cluster <- fct_recode(
    sce_transcript_marker$cluster,
    `cyc_rad_glia.1` = "1",
    `gaba.2` = "2",
    `gaba.3` = "3",
    `glut.4` = "4",
    `glut.5` = "5",
    `rad_glia.6` = "6",
    `glut.7` = "7",
    `glut.8` = "8",
    `glut.9` = "9",
    `cr.10` = "10",
    `glut.11` = "11"
)

cluster_comparisons <- combn(levels(sce_transcript_marker$cluster), 2,
                             simplify = FALSE)
marker_results_list <- bplapply(cluster_comparisons, function(compared_clusters) {
    findMarkers(sce_transcript_marker,
                groups = sce_transcript_marker$cluster,
                restrict = c(compared_clusters[1], compared_clusters[2]),
                test.type = "wilcox",
                pval.type = "any",
                direction = "up",
                row.data = rowData(sce_transcript_marker))
}, BPPARAM = BPPARAM)
names(marker_results_list) <- sapply(cluster_comparisons, paste0, collapse = "__")

marker_df_list <- lapply(names(marker_results_list), function(contrast_name) {
    marker_results <- marker_results_list[[contrast_name]]
    cluster_1_df <- marker_results[[1]] %>%
        as.data.frame() %>%
        transmute(
            transcript_id = transcript_id,
            compatible_tx = compatible_tx,
            gene_id = gene_id,
            gene_name = gene_name,
            pvalue = p.value,
            fdr = FDR,
            auc = summary.AUC,
            cluster_1 = names(marker_results)[1],
            cluster_2 = names(marker_results)[2],
            contrast = contrast_name
        ) %>%
        `rownames<-`(NULL)
    cluster_2_df <- marker_results[[2]] %>%
        as.data.frame() %>%
        transmute(
            transcript_id = transcript_id,
            compatible_tx = compatible_tx,
            gene_id = gene_id,
            gene_name = gene_name,
            pvalue = p.value,
            fdr = FDR,
            auc = summary.AUC,
            cluster_1 = names(marker_results)[2],
            cluster_2 = names(marker_results)[1],
            contrast = contrast_name
        ) %>%
        `rownames<-`(NULL)
    return(rbind(cluster_1_df, cluster_2_df))
})
marker_df <- do.call(rbind, marker_df_list)
saveRDS(marker_df, file.path(result_dir, "marker_df.rds"))

################################################################################

# Detect isoform swicthing events

isoswitch_df <- marker_df %>%
    filter(fdr <= 0.05) %>%
    group_by(gene_id, contrast) %>%
    mutate(
        n_clusters = length(unique(cluster_1)),
        n_transcripts = length(unique(transcript_id))
    ) %>%
    ungroup() %>%
    filter(n_clusters > 1, n_transcripts > 1) %>%
    select(-n_clusters, -n_transcripts) %>%
    arrange(contrast, gene_id)
saveRDS(isoswitch_df, file.path(result_dir, "isoswitch_df.rds"))

write.csv(
    isoswitch_df,
    "results/isoform_switch.csv",
    row.names = FALSE
)

################################################################################

# Create gene & transcript expression UMAP plots for significant results

isoswitch_data <- isoswitch_df[, 1:4] %>%
    distinct() %>%
    arrange(gene_name)

plot_isoform_umap_expression <- function(transcript_id, compatible_tx,
                                         gene_id, gene_name) {
    p1 <- plot_density(sce, transcript_id,
                       size = 1.5) +
        labs(title = "Transcript expression density") +
        theme(legend.position = "right",
              plot.title = element_text(size = 13),
              legend.title = element_text(size = 11))
    p2 <- plot_density(sce, gene_id,
                       size = 1.5) +
        labs(title = "Gene expression density") +
        theme(legend.position = "right",
              plot.title = element_text(size = 13),
              legend.title = element_text(size = 11))
    p3 <- dittoDimPlot(sce,
                       transcript_id,
                       reduction.use = "UMAP",
                       size = 1.5, main = "Transcript expression",
                       legend.title = "Logcounts",
                       order = "increasing")
    p4 <- dittoDimPlot(sce, gene_id,
                       reduction.use = "UMAP",
                       size = 1.5, main = "Gene expression",
                       legend.title = "Logcounts",
                       order = "increasing")

    patchwork <- (p1 + p2) / (p3 + p4)
    patchwork + plot_annotation(title = glue("{transcript_id} ({compatible_tx})"),
                                subtitle = glue("{gene_name} ({gene_id})")) &
        theme(plot.title = element_text(size = 12))
}

dir.create("results/temp", recursive = TRUE)
pdf_files <- sapply(seq(nrow(isoswitch_data)), function(i) {
    umap_plot <- plot_isoform_umap_expression(isoswitch_data$transcript_id[i],
                                              isoswitch_data$compatible_tx[i],
                                              isoswitch_data$gene_id[i],
                                              isoswitch_data$gene_name[i])
    pdf_file <- glue("results/temp/{i}.pdf")
    ggsave(pdf_file, umap_plot, height = 5, width = 8)
    return(pdf_file)
})
pdf_combine(pdf_files, output = "results/umap_isoswitch.pdf")
unlink("results/temp", recursive = TRUE)

################################################################################

# Print session information
sessionInfo()
