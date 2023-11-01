#! /usr/bin/env bash

# Run the R scripts
singularity exec ../singularity/isosceles.sif Rscript 00_run_isosceles.R
singularity exec ../singularity/isosceles.sif Rscript 01_scrnaseq_analysis.R
singularity exec ../singularity/isosceles.sif Rscript 02_run_dexseq_intra.R
singularity exec ../singularity/isosceles.sif Rscript 03_run_dexseq_inter.R
singularity exec ../singularity/isosceles.sif Rscript 04_psi_heatmap_intra.R
singularity exec ../singularity/isosceles.sif Rscript 05_psi_heatmap_inter.R
singularity exec ../singularity/isosceles.sif Rscript 06_other_plots.R
singularity exec ../singularity/isosceles.sif Rscript 07_isoswitch_analysis.R
