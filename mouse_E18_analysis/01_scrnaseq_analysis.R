#!/usr/bin/env Rscript

suppressMessages({
    library(tidyverse)
    library(glue)
    library(Isosceles)
    library(scran)
    library(scater)
    library(igraph)
    library(bluster)
    library(slingshot)
    library(BiocParallel)
    library(BiocNeighbors)
})

# Set the number of CPUs/threads for the analysis
ncpu <- 1

# Global parameters
BNPARAM <- KmknnParam()
BPPARAM <- MulticoreParam(ncpu)

result_dir <- "01_scrnaseq_analysis"
dir.create(result_dir, recursive = TRUE)

# Read Isosceles SE objects
se_tcc <- readRDS("00_run_isosceles/se_tcc.rds")
se_gene <- readRDS("00_run_isosceles/se_gene.rds")
se_transcript <- readRDS("00_run_isosceles/se_transcript.rds")
se_psi <- readRDS("00_run_isosceles/se_psi.rds")

################################################################################

# Filtering the data

# Keep only cells with at least 100 expressed genes
cell_selector <- colSums(assay(se_gene, "counts") > 0) >= 100
se_tcc <- se_tcc[, cell_selector]
se_gene <- se_gene[, cell_selector]
se_transcript <- se_transcript[, cell_selector]
se_psi <- se_psi[, cell_selector]

# Assign batch labels to cells and keep only the 951 cells sequencing data
cell_labels <- sapply(strsplit(colnames(se_gene), "\\."), "[", 1)
se_tcc$label <- cell_labels
se_gene$label <- cell_labels
se_transcript$label <- cell_labels
se_psi$label <- cell_labels
cell_selector <- se_gene$label == "951_cells"
se_tcc <- se_tcc[, cell_selector]
se_gene <- se_gene[, cell_selector]
se_transcript <- se_transcript[, cell_selector]
se_psi <- se_psi[, cell_selector]

################################################################################

# Basic scRNA-Seq analysis

# Prepare the joint gene & transcript SCE object
sce <- SingleCellExperiment(
    assays = list(counts = rbind(
        assay(se_gene, "counts"),
        assay(se_transcript, "counts")
    ))
)

# Detecting highly variable features
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)
dec <- modelGeneVar(sce)
top_hvgs <- getTopHVGs(dec, n = 4000)
saveRDS(top_hvgs, file.path(result_dir, "top_hvgs.rds"))
# 3,760 features (1,735 genes and 2,025 transcripts)
# Identified features include every result with dec$bio > 0

# Dimensionality reduction
set.seed(42)
sce <- runPCA(sce, ncomponents = 30, subset_row = top_hvgs, scale = TRUE)
set.seed(42)
sce <- runUMAP(sce, dimred = "PCA")

# Cell clustering
set.seed(42)
snn_graph <- buildSNNGraph(sce, use.dimred = "PCA", k = 10,
                           BPPARAM = BPPARAM, BNPARAM = BNPARAM)
set.seed(42)
sce$cluster <- as.factor(cluster_louvain(
    snn_graph, resolution = 2
)$membership)
saveRDS(sce, file.path(result_dir, "sce.rds"))

# Process the other SCE objects
sce_gene <- as(se_gene, "SingleCellExperiment")
sce_gene <- computeLibraryFactors(sce_gene)
sce_gene <- logNormCounts(sce_gene)
reducedDim(sce_gene, "UMAP") <- reducedDim(sce, "UMAP")
sce_gene$cluster <- sce$cluster
saveRDS(sce_gene, file.path(result_dir, "sce_gene.rds"))

sce_transcript <- as(se_transcript, "SingleCellExperiment")
sce_transcript <- computeLibraryFactors(sce_transcript)
sce_transcript <- logNormCounts(sce_transcript)
reducedDim(sce_transcript, "UMAP") <- reducedDim(sce_gene, "UMAP")
sce_transcript$cluster <- sce$cluster
saveRDS(sce_transcript, file.path(result_dir, "sce_transcript.rds"))

sce_psi <- as(se_psi, "SingleCellExperiment")
reducedDim(sce_psi, "UMAP") <- reducedDim(sce_gene, "UMAP")
sce_psi$cluster <- sce$cluster
saveRDS(sce_psi, file.path(result_dir, "sce_psi.rds"))

################################################################################

# Preparing pseudobulk SE objects

se_pseudobulk_tcc <- pseudobulk_tcc(
    se_tcc = se_tcc,
    cell_labels = sce$cluster
)
saveRDS(se_pseudobulk_tcc, file.path(result_dir, "se_pseudobulk_tcc.rds"))

se_pseudobulk_gene <- tcc_to_gene(
    se_tcc = se_pseudobulk_tcc
)
saveRDS(se_pseudobulk_gene, file.path(result_dir, "se_pseudobulk_gene.rds"))

se_pseudobulk_transcript <- tcc_to_transcript(
    se_tcc = se_pseudobulk_tcc, em.maxiter = 250, em.conv = 0.01,
    use_length_normalization = FALSE, ncpu = ncpu
)
saveRDS(se_pseudobulk_transcript, file.path(result_dir, "se_pseudobulk_transcript.rds"))

se_pseudobulk_psi <- transcript_to_psi(se_pseudobulk_transcript, ncpu = ncpu)
saveRDS(se_pseudobulk_psi, file.path(result_dir, "se_pseudobulk_psi.rds"))

################################################################################

# Pseudotime analysis using Slingshot

# Assign clusters to cell types / trajectories
clusters_glut <- c("4", "5", "7", "8",  "9", "11")
clusters_gaba <- c("2", "3")
clusters_rad_glia <- "6"
clusters_cyc_rad_glia <- "1"
clusters_cr <- "10"

# Glutamatergic neurons
sce$cluster_glut <- ifelse(sce$cluster %in% clusters_glut,
                           sce$cluster, "-1")
set.seed(42)
sce_glut <- slingshot(sce, reducedDim = "PCA", use.median = TRUE,
                      clusterLabels = sce$cluster_glut,
                      start.clus = "4",
                      stretch = 0, approx_points = 200)

# GABAergic neurons
sce$cluster_gaba <- ifelse(sce$cluster %in% clusters_gaba,
                           sce$cluster, "-1")
set.seed(42)
sce_gaba <- slingshot(sce, reducedDim = "PCA",
                      clusterLabels = sce$cluster_gaba,
                      start.clus = "2",
                      stretch = 0, approx_points = 200)

# Radial glia
sce$cluster_rad_glia <- ifelse(sce$cluster %in% clusters_rad_glia,
                               sce$cluster, "-1")
set.seed(42)
sce_rad_glia <- slingshot(sce, reducedDim = "PCA", use.median = TRUE,
                          clusterLabels = sce$cluster_rad_glia,
                          start.clus = "6",
                          stretch = 0, approx_points = 200)

# Cycling radial glia
sce$cluster_cyc_rad_glia <- ifelse(sce$cluster %in% clusters_cyc_rad_glia,
                                   sce$cluster, "-1")
set.seed(42)
sce_cyc_rad_glia <- slingshot(sce, reducedDim = "PCA", use.median = TRUE,
                              clusterLabels = sce$cluster_cyc_rad_glia,
                              start.clus = "1",
                              stretch = 0, approx_points = 200)

# Cajal-Retzius cells
sce$cluster_cr <- ifelse(sce$cluster %in% clusters_cr,
                         sce$cluster, "-1")
set.seed(42)
sce_cr <- slingshot(sce, reducedDim = "PCA", use.median = TRUE,
                    clusterLabels = sce$cluster_cr,
                    start.clus = "10",
                    stretch = 0, approx_points = 200)

# Prepare the pseudotime matrix
pseudotime_matrix <- cbind(
    glut_1 = colData(sce_glut)$slingPseudotime_1,
    glut_2 = colData(sce_glut)$slingPseudotime_2,
    gaba = colData(sce_gaba)$slingPseudotime_1,
    rad_glia = colData(sce_rad_glia)$slingPseudotime_1,
    cyc_rad_glia = colData(sce_cyc_rad_glia)$slingPseudotime_1,
    cr = colData(sce_cr)$slingPseudotime_1
)
rownames(pseudotime_matrix) <- colnames(sce_glut)
saveRDS(pseudotime_matrix, file.path(result_dir, "pseudotime_matrix.rds"))

################################################################################

# Prepare PSI regions of interest for each trajectory
psi_events_list <- lapply(colnames(pseudotime_matrix), function(traj_name) {
    pseudotime <- pseudotime_matrix[, traj_name]
    pseudotime <- pseudotime[!is.na(pseudotime)]
    se_psi_traj <- se_psi[, names(pseudotime)]
    psi_assay <- assay(se_psi_traj, "psi")
    psi_assay <- psi_assay[rowData(se_psi_traj)$gene_id %in% top_hvgs,]
    feature_selector <- (!grepl(":TSS", rownames(psi_assay))) &
        (!grepl(":TES", rownames(psi_assay)))
    psi_assay <- psi_assay[feature_selector,]
    psi_selector <- (apply(psi_assay, 1, mean) > 0.025) &
        (apply(psi_assay, 1, mean) < 0.975) &
        (apply(psi_assay, 1, function(x) { sum((x > 0) & (x < 0.999) & (x != 0.5)) }) > 30) &
        (apply(psi_assay, 1, function(x) { sum(x > 0.1) }) > 30)
    psi_assay <- psi_assay[psi_selector,]
    psi_events <- rownames(psi_assay)
    return(psi_events)
})
names(psi_events_list) <- colnames(pseudotime_matrix)
saveRDS(psi_events_list, file.path(result_dir, "psi_events_list.rds"))

################################################################################

# Print session information
sessionInfo()
