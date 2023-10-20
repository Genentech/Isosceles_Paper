#!/usr/bin/env Rscript

suppressMessages({
    library(tidyverse)
    library(scales)
    library(glue)
    library(ggnewscale)
    library(scran)
    library(scater)
    library(BiocParallel)
    library(RColorBrewer)
    library(pheatmap)
    library(dittoSeq)
    library(Nebulosa)
})

# Set the number of CPUs/threads for the analysis
ncpu <- 1

# Global parameters
BPPARAM <- MulticoreParam(ncpu)

result_dir <- "06_other_plots"
dir.create(result_dir, recursive = TRUE)

# Read the analysis results from previous steps
sce <- readRDS("01_scrnaseq_analysis/sce.rds")
sce_gene <- readRDS("01_scrnaseq_analysis/sce_gene.rds")
se_pseudobulk_gene <- readRDS("01_scrnaseq_analysis/se_pseudobulk_gene.rds")
pseudotime_matrix <- readRDS("01_scrnaseq_analysis/pseudotime_matrix.rds")
heatmap_window_pseudotime_list <- readRDS("04_psi_heatmap_intra/heatmap_window_pseudotime_list.rds")

################################################################################

# Create the UMAP cluster plot

sce_plot <- sce
sce_plot$cluster <- fct_recode(
    sce_plot$cluster,
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

umap_plot <- dittoDimPlot(sce_plot, "cluster",
                          reduction.use = "UMAP", main = "",
                          size = 2, do.label = TRUE, labels.size = 4,
                          labels.highlight = TRUE, legend.show = FALSE) +
    scale_colour_manual(
        values = c(
            glut.4 = brewer.pal(n = 7, name = "YlGnBu")[3],
            glut.8 = brewer.pal(n = 7, name = "YlGnBu")[4],
            glut.5 = brewer.pal(n = 7, name = "YlGnBu")[5],
            glut.7 = brewer.pal(n = 7, name = "YlGnBu")[6],
            glut.11 = brewer.pal(n = 7, name = "YlGn")[5],
            glut.9 = brewer.pal(n = 7, name = "YlGn")[6],
            gaba.2 = brewer.pal(n = 7, name = "Reds")[4],
            gaba.3 = brewer.pal(n = 7, name = "Reds")[6],
            rad_glia.6 = colorRampPalette(c("white", "deeppink"))(7)[6],
            cyc_rad_glia.1 = colorRampPalette(c("white", "purple"))(7)[6],
            cr.10 = colorRampPalette(c("white", "gold"))(7)[6]
        )
    )
ggsave("results/umap_cluster.pdf", umap_plot, height = 5, width = 7)

################################################################################

# Create the UMAP trajectory plot

traj_colors <- list(
    glut_1 = brewer.pal(n = 7, name = "YlGnBu")[2:6],
    glut_2 = brewer.pal(n = 7, name = "YlGn")[2:6],
    gaba = brewer.pal(n = 7, name = "Reds")[2:6],
    rad_glia = c("white", "deeppink"),
    cyc_rad_glia = c("white", "purple"),
    cr = c("white", "gold")
)

pseudotime_lim <- sapply(heatmap_window_pseudotime_list, function(x) {range(x)})
rownames(pseudotime_lim) <- c("min", "max")

umap_df <- as.data.frame(reducedDim(sce, "UMAP"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    pseudotime <- pseudotime_matrix[, traj_name]
    pseudotime <- pseudotime[!is.na(pseudotime)]
    umap_traj_df <- umap_df[names(pseudotime),]
    umap_traj_df$pseudotime <- pseudotime
    return(umap_traj_df)
})
names(umap_df_list) <- colnames(pseudotime_matrix)

cluster_umap_df <- umap_df %>%
    mutate(cluster = sce$cluster) %>%
    group_by(cluster) %>%
    summarise(
        UMAP_1 = mean(UMAP_1),
        UMAP_2 = mean(UMAP_2)
    )

umap_lineages <- list(
    glut_1 = c("4", "8", "5", "7"),
    glut_2 = c("4", "11", "9"),
    gaba = c("2", "3")
)

lineage_umap_df <- umap_lineages %>%
    enframe() %>%
    unchop(value) %>%
    dplyr::rename(lineage = name, cluster = value) %>%
    left_join(cluster_umap_df)

umap_plot <- ggplot() +
    geom_point(
        data = umap_df_list[["glut_1"]],
        shape = 21, colour = "grey", size = 3, stroke = 0,
        mapping = aes(x = UMAP_1, y = UMAP_2, fill = pseudotime)
    ) +
    scale_fill_gradientn(
        colours = traj_colors[["glut_1"]],
        limits = pseudotime_lim[, "glut_1"],
        oob = scales::squish
    ) +
    new_scale_fill() +
    geom_point(
        data = umap_df_list[["glut_2"]],
        shape = 21, colour = "grey", size = 3, stroke = 0,
        mapping = aes(x = UMAP_1, y = UMAP_2, fill = pseudotime)
    ) +
    scale_fill_gradientn(
        colours = traj_colors[["glut_2"]],
        limits = pseudotime_lim[, "glut_2"],
        oob = scales::squish
    ) +
    new_scale_fill() +
    geom_point(
        data = umap_df_list[["gaba"]],
        shape = 21, colour = "grey", size = 3, stroke = 0,
        mapping = aes(x = UMAP_1, y = UMAP_2, fill = pseudotime)
    ) +
    scale_fill_gradientn(
        colours = traj_colors[["gaba"]],
        limits = pseudotime_lim[, "gaba"],
        oob = scales::squish
    ) +
    new_scale_fill() +
    geom_point(
        data = umap_df_list[["rad_glia"]],
        shape = 21, colour = "grey", size = 3, stroke = 0,
        mapping = aes(x = UMAP_1, y = UMAP_2, fill = pseudotime)
    ) +
    scale_fill_gradientn(
        colours = traj_colors[["rad_glia"]],
        limits = pseudotime_lim[, "rad_glia"],
        oob = scales::squish
    ) +
    new_scale_fill() +
    geom_point(
        data = umap_df_list[["cyc_rad_glia"]],
        shape = 21, colour = "grey", size = 3, stroke = 0,
        mapping = aes(x = UMAP_1, y = UMAP_2, fill = pseudotime)
    ) +
    scale_fill_gradientn(
        colours = traj_colors[["cyc_rad_glia"]],
        limits = pseudotime_lim[, "cyc_rad_glia"],
        oob = scales::squish
    ) +
    new_scale_fill() +
    geom_point(
        data = umap_df_list[["cr"]],
        shape = 21, colour = "grey", size = 3, stroke = 0,
        mapping = aes(x = UMAP_1, y = UMAP_2, fill = pseudotime)
    ) +
    scale_fill_gradientn(
        colours = traj_colors[["cr"]],
        limits = pseudotime_lim[, "cr"],
        oob = scales::squish
    ) +
    geom_path(
        data = lineage_umap_df,
        size = 0.8,
        mapping = aes(x = UMAP_1, y = UMAP_2, group = lineage),
        arrow = arrow(ends = "last",
                      type = "closed",
                      length = unit(0.1, "inches"))
    ) +
    theme_classic() +
    theme(legend.position = "none")

ggsave("results/umap_trajectory.pdf", umap_plot, height = 5, width = 7)

################################################################################

# Create the marker gene heatmap

traj_colors <- c(
    glut_1.4 = brewer.pal(n = 7, name = "YlGnBu")[3],
    glut_1.8 = brewer.pal(n = 7, name = "YlGnBu")[4],
    glut_1.5 = brewer.pal(n = 7, name = "YlGnBu")[5],
    glut_1.7 = brewer.pal(n = 7, name = "YlGnBu")[6],
    glut_2.4 = brewer.pal(n = 7, name = "YlGn")[4],
    glut_2.11 = brewer.pal(n = 7, name = "YlGn")[5],
    glut_2.9 = brewer.pal(n = 7, name = "YlGn")[6],
    gaba.2 = brewer.pal(n = 7, name = "Reds")[4],
    gaba.3 = brewer.pal(n = 7, name = "Reds")[6],
    rad_glia.6 = colorRampPalette(c("white", "deeppink"))(7)[6],
    cyc_rad_glia.1 = colorRampPalette(c("white", "purple"))(7)[6],
    cr.10 = colorRampPalette(c("white", "gold"))(7)[6]
)

column_order <- c(4, 8, 5, 7, 4, 11, 9, 2, 3, 6, 1, 10)
traj_names <- rep(
    c("glut_1", "glut_2", "gaba", "rad_glia", "cyc_rad_glia", "cr"),
    c(4, 3, 2, 1, 1, 1)
)
mat_colnames <- paste0(traj_names, ".", column_order)
col_ann_df <- data.frame(
    trajectory = mat_colnames # traj_names
)
rownames(col_ann_df) <- mat_colnames

sce_gene_marker <- sce_gene
rownames(sce_gene_marker) <- make.unique(
    rowData(sce_gene_marker)$gene_name, sep = "@"
)
se_pseudobulk_gene_marker <- se_pseudobulk_gene
rownames(se_pseudobulk_gene_marker) <- make.unique(
    rowData(se_pseudobulk_gene_marker)$gene_name, sep = "@"
)

marker_genes <- findMarkers(sce_gene_marker,
                            groups = sce_gene_marker$cluster,
                            test.type = "t",
                            pval.type = "some",
                            direction = "up",
                            BPPARAM = BPPARAM)
saveRDS(marker_genes, file.path(result_dir, "marker_genes.rds"))

top_marker_genes <- lapply(
        marker_genes, function(x) rownames(x[x$FDR <= 0.05,])[1:5]
    ) %>%
    unlist() %>%
    unique() %>%
    sort()

mat_gene <- log2(assay(se_pseudobulk_gene_marker, "tpm") + 1)
mat_gene <- mat_gene[top_marker_genes,]
mat_gene <- mat_gene[, column_order]
colnames(mat_gene) <- mat_colnames

heatmap_gene <- pheatmap(
    mat_gene,
    cluster_rows = TRUE, cluster_cols = FALSE,
    show_rownames = TRUE, show_colnames = FALSE,
    annotation_col = col_ann_df,
    annotation_colors = list(trajectory = traj_colors),
    annotation_legend = FALSE, annotation_names_col = FALSE,
    fontsize_row = 6,
    scale = "row"
)
ggsave("results/heatmap_marker_genes.pdf", heatmap_gene,
       height = 7, width = 7)

################################################################################

# Create the marker gene set score plots

traj_colors <- c(
    glut_1.4 = brewer.pal(n = 7, name = "YlGnBu")[3],
    glut_1.8 = brewer.pal(n = 7, name = "YlGnBu")[4],
    glut_1.5 = brewer.pal(n = 7, name = "YlGnBu")[5],
    glut_1.7 = brewer.pal(n = 7, name = "YlGnBu")[6],
    glut_2.4 = brewer.pal(n = 7, name = "YlGn")[4],
    glut_2.11 = brewer.pal(n = 7, name = "YlGn")[5],
    glut_2.9 = brewer.pal(n = 7, name = "YlGn")[6],
    gaba.2 = brewer.pal(n = 7, name = "Reds")[4],
    gaba.3 = brewer.pal(n = 7, name = "Reds")[6],
    rad_glia.6 = colorRampPalette(c("white", "deeppink"))(7)[6],
    cyc_rad_glia.1 = colorRampPalette(c("white", "purple"))(7)[6],
    cr.10 = colorRampPalette(c("white", "gold"))(7)[6]
)

column_order <- c(4, 8, 5, 7, 4, 11, 9, 2, 3, 6, 1, 10)
traj_names <- rep(
    c("glut_1", "glut_2", "gaba", "rad_glia", "cyc_rad_glia", "cr"),
    c(4, 3, 2, 1, 1, 1)
)
mat_colnames <- paste0(traj_names, ".", column_order)
col_ann_df <- data.frame(
    trajectory = mat_colnames # traj_names
)
rownames(col_ann_df) <- mat_colnames

marker_sets <- list(
    Progenitors = c("Neurog2", "Eomes"),
    CR_cells = c("Snhg11", "Lhx5", "Reln"),
    GABA = c("Gad2", "Gad1", "Dlx2"),
    Mat_GABA = c("Maf", "Mafb", "Arx"),
    Glut = c("Neurod6", "Neurod2"),
    Mat_Glut = c("Camk2b", "Opcml", "Crym"),
    Rad_glia = c("Fabp7", "Vim", "Dbi"),
    Cyc_rad_glia = c("Cenpf", "Top2a", "Mki67")
)

se_pseudobulk_gene_marker <- se_pseudobulk_gene
rownames(se_pseudobulk_gene_marker) <- make.unique(
    rowData(se_pseudobulk_gene_marker)$gene_name, sep = "@"
)

marker_scores <- sumCountsAcrossFeatures(se_pseudobulk_gene_marker,
                                         marker_sets,
                                         exprs_values = "tpm",
                                         average = TRUE)
saveRDS(marker_scores, file.path(result_dir, "marker_scores.rds"))

mat_gene_set <- log2(marker_scores + 1)
mat_gene_set <- mat_gene_set[, column_order]
colnames(mat_gene_set) <- mat_colnames

heatmap_gene_set <- pheatmap(
    mat_gene_set,
    cluster_rows = TRUE, cluster_cols = FALSE,
    show_rownames = TRUE, show_colnames = FALSE,
    annotation_col = col_ann_df,
    annotation_colors = list(trajectory = traj_colors),
    annotation_legend = FALSE, annotation_names_col = FALSE,
    fontsize_row = 8,
    scale = "row"
)
ggsave("results/heatmap_marker_gene_sets_000.pdf", heatmap_gene_set,
       height = 5, width = 7)

sce_gene_marker <- sce_gene
rownames(sce_gene_marker) <- make.unique(
    rowData(sce_gene_marker)$gene_name, sep = "@"
)

cell_marker_scores <- sumCountsAcrossFeatures(sce_gene_marker,
                                              marker_sets,
                                              exprs_values = "logcounts",
                                              average = TRUE)

sce_gene_set <- SingleCellExperiment(
    assays = list(logcounts = cell_marker_scores),
    colData =  colData(sce_gene_marker)
)
reducedDim(sce_gene_set, "UMAP") <- reducedDim(sce_gene_marker, "UMAP")
saveRDS(sce_gene_set, file.path(result_dir, "sce_gene_set.rds"))

density_plot <- plot_density(
    sce_gene_set,
    c("Progenitors", "Glut", "Mat_Glut", "GABA", "Mat_GABA",
      "Rad_glia", "Cyc_rad_glia", "CR_cells"),
    slot = "logcounts", size = 1
)
ggsave("results/density_plot_marker_gene_sets.pdf", density_plot,
       height = 7, width = 10)

################################################################################

# Print session information
sessionInfo()
