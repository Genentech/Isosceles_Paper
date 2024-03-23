# Preparing the benchmark reports

## Content

  * **run.sh** - the script for building the reports

## Report list

  * **simulated_bulk_benchmarks.ipynb** - simulated bulk RNA-Seq data benchmarks
  * **simulated_sc_benchmarks.Rmd** - simulated scRNA-Seq data benchmarks
  * **nanopore_bulk_sc_benchmarks_hvts.Rmd** - Nanopore bulk RNA-Seq and scRNA-Seq data benchmarks (top 4000 highly variable transcripts)
  * **nanopore_bulk_sc_benchmarks_tpm.Rmd** - Nanopore bulk RNA-Seq and scRNA-Seq data benchmarks (all >= 1 TPM transcripts)
  * **nanopore_bulk_igrov_mean_rel_diff.Rmd** - Nanopore IGROV-1 bulk RNA-Seq mean relative difference benchmarks
  * **nanopore_bulk_sc_igrov_mean_rel_diff.Rmd** - Nanopore IGROV-1 bulk RNA-Seq and scRNA-Seq mean relative difference benchmarks
  * **nanopore_bulk_sc_igrov_cor.Rmd** - Nanopore IGROV-1 bulk RNA-Seq and scRNA-Seq correlation benchmarks
  * **nanopore_sc_umap.Rmd** - Nanopore scRNA-Seq UMAP plots
  * **sirv_benchmarks.Rmd** - SIRV data benchmarks
  * **sequin_benchmarks.Rmd** - Sequin data benchmarks

## Required programs and libraries

  * Singularity (>= 3.5.1)
