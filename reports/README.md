# Preparing the benchmark reports

## Content

  * **run.sh** - the script for building the reports

## Report list

  * **simulated_bulk_benchmarks.ipynb** - simulated bulk RNA-Seq data benchmarks
  * **simulated_sc_benchmarks.Rmd** - simulated scRNA-Seq data benchmarks
  * **nanopore_bulk_sc_benchmarks_hvts_\*.Rmd** - Nanopore bulk RNA-Seq and scRNA-Seq data benchmarks for different numbers of top highly variable transcripts
  * **nanopore_bulk_sc_benchmarks_hvts_data/** - a directory for storing Nanopore bulk RNA-Seq and scRNA-Seq data benchmarks results
  * **nanopore_bulk_sc_benchmarks_hvts_summary.Rmd** - Nanopore bulk RNA-Seq and scRNA-Seq data benchmarks (top highly variable transcripts results summary)
  * **nanopore_bulk_sc_benchmarks_tpm.Rmd** - Nanopore bulk RNA-Seq and scRNA-Seq data benchmarks (all >= 1 TPM transcripts)
  * **simulated_ovarian_bulk_sc_benchmarks_hvts_\*.Rmd** - Simulated ovarian cell line bulk RNA-Seq and scRNA-Seq data benchmarks for different numbers of top highly variable transcripts
  * **simulated_ovarian_bulk_sc_benchmarks_hvts_data/** - a directory for storing simulated ovarian cell line bulk RNA-Seq and scRNA-Seq data benchmarks results
  * **simulated_ovarian_bulk_sc_benchmarks_hvts_summary.Rmd** - Simulated ovarian cell line bulk RNA-Seq and scRNA-Seq data benchmarks (top highly variable transcripts results summary)
  * **simulated_ovarian_bulk_sc_benchmarks_tpm.Rmd** - Simulated ovarian cell line bulk RNA-Seq and scRNA-Seq data benchmarks (all >= 1 TPM transcripts)
  * **nanopore_sc_umi_benchmarks.Rmd** - Nanopore scRNA-Seq Sicelore UMI vs wf-single-cell UMI benchmarks
  * **nanopore_bulk_sc_umi_benchmarks_hvts_4000.Rmd** - Nanopore bulk RNA-Seq and scRNA-Seq Sicelore UMI vs wf-single-cell UMI benchmarks (4000 top highly variable transcripts)
  * **nanopore_bulk_sc_umi_benchmarks_tpm.Rmd** - Nanopore bulk RNA-Seq and scRNA-Seq Sicelore UMI vs wf-single-cell UMI benchmarks (all >= 1 TPM transcripts)
  * **nanopore_bulk_igrov_mean_rel_diff.Rmd** - Nanopore IGROV-1 bulk RNA-Seq mean relative difference benchmarks
  * **nanopore_bulk_sc_igrov_mean_rel_diff.Rmd** - Nanopore IGROV-1 bulk RNA-Seq and scRNA-Seq mean relative difference benchmarks
  * **nanopore_bulk_sc_igrov_cor.Rmd** - Nanopore IGROV-1 bulk RNA-Seq and scRNA-Seq correlation benchmarks
  * **nanopore_sc_umap.Rmd** - Nanopore scRNA-Seq UMAP plots
  * **sirv_benchmarks.Rmd** - SIRV data benchmarks
  * **sequin_benchmarks.Rmd** - Sequin data benchmarks
  * **pacbio_benchmarks.Rmd** - Pacbio, Nanopore and Illumina GM12878 data benchmarks

## Required programs and libraries

  * Singularity (>= 3.5.1)
