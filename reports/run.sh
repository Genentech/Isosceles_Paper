#! /usr/bin/env bash

# Build the Jupyter reports
singularity exec ../singularity/reports_python.sif \
  jupyter nbconvert --execute --ExecutePreprocessor.timeout=-1 \
  --to notebook --inplace simulated_bulk_benchmarks.ipynb

# Build the R reports
singularity exec ../singularity/reports_r.sif Rscript -e 'rmarkdown::render("simulated_sc_benchmarks.Rmd")'
singularity exec ../singularity/reports_r.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_benchmarks_hvts.Rmd")'
singularity exec ../singularity/reports_r.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_benchmarks_tpm.Rmd")'
singularity exec ../singularity/reports_r.sif Rscript -e 'rmarkdown::render("nanopore_bulk_igrov_mean_rel_diff.Rmd")'
singularity exec ../singularity/reports_r.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_igrov_mean_rel_diff.Rmd")'
singularity exec ../singularity/reports_r.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_igrov_cor.Rmd")'
singularity exec ../singularity/reports_r.sif Rscript -e 'rmarkdown::render("nanopore_sc_umap.Rmd")'
