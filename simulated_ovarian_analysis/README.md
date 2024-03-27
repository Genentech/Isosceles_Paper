# Simulated ovarian cell line RNA-Seq data analysis

## Content

  * **run.sh** - the script for running the analysis
  * **run_isosceles_bulk.R** - an R script used for running Isosceles (bulk RNA-Seq data analysis)
  * **run_isosceles_sc.R** - an R script used for running Isosceles (scRNA-Seq data analysis)
  * **add_bam_tags.py** - a Python script used for adding scRNA-Seq tags to the bulk RNA-Seq BAM files
  * **prepare_report_data.sh** - the script for preparing report data from run.sh script results
  * **download_report_data.sh** - the script for downloading the already prepared report data

## Required programs and libraries

  * conda
  * Singularity (>= 3.5.1)
  * samtools (>= 1.10)
  * bedtools
  * UCSC Genome Browser utilities
  * Java environment (>= 1.8)
