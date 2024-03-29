# Nanopore scRNA-Seq data analysis

## Content

  * **run.sh** - the script for running the analysis
  * **run_isosceles.R** - an R script used for running Isosceles (Sicelore input)
  * **run_isosceles_wf.R** - an R script used for running Isosceles (wf-single-cell input)
  * **prepare_report_data.sh** - the script for preparing report data from run.sh script results
  * **download_report_data.sh** - the script for downloading the already prepared report data

## Required programs and libraries

  * conda
  * Singularity (>= 3.5.1)
  * samtools (>= 1.18)
  * UCSC Genome Browser utilities
  * Java environment (>= 1.8)
  * Cell Ranger (>= 7.1.0)
  * Nextflow (>= 23.04.4)
