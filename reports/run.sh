#! /usr/bin/env bash

# Build the Jupyter reports
singularity exec ../singularity/reports_python.sif \
  jupyter nbconvert --execute --ExecutePreprocessor.timeout=-1 \
  --to notebook --inplace simulated_bulk_benchmarks.ipynb

# Build the R reports
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("simulated_sc_benchmarks.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_benchmarks_hvts_500.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_benchmarks_hvts_1000.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_benchmarks_hvts_2000.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_benchmarks_hvts_4000.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_benchmarks_hvts_6000.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_benchmarks_hvts_10000.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_benchmarks_hvts_summary.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_benchmarks_tpm.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("simulated_ovarian_bulk_sc_benchmarks_hvts_500.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("simulated_ovarian_bulk_sc_benchmarks_hvts_1000.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("simulated_ovarian_bulk_sc_benchmarks_hvts_2000.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("simulated_ovarian_bulk_sc_benchmarks_hvts_4000.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("simulated_ovarian_bulk_sc_benchmarks_hvts_6000.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("simulated_ovarian_bulk_sc_benchmarks_hvts_10000.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("simulated_ovarian_bulk_sc_benchmarks_hvts_summary.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("simulated_ovarian_bulk_sc_benchmarks_tpm.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_sc_umi_benchmarks.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_umi_benchmarks_tpm.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_umi_benchmarks_hvts_4000.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_igrov_mean_rel_diff.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_igrov_mean_rel_diff.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_bulk_sc_igrov_cor.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("nanopore_sc_umap.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("sirv_benchmarks.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("sequin_benchmarks.Rmd")'
singularity exec ../singularity/isosceles.sif Rscript -e 'rmarkdown::render("pacbio_benchmarks.Rmd")'
