# Nanopore bulk RNA-Seq data analysis

## Content

  * **run.sh** - the script for running the analysis
  * **run_isosceles.R** - an R script used for running Isosceles
  * **run_bambu.R** - an R script used for running bambu
  * **add_bam_tags.py** - a Python script used for adding scRNA-Seq tags to the bulk RNA-Seq BAM files
  * **prepare_report_data.sh** - the script for preparing report data from run.sh script results
  * **download_report_data.sh** - the script for downloading the already prepared report data

## Required programs and libraries

  * conda
  * Singularity (>= 3.5.1)
  * samtools (>= 1.10)
  * UCSC Genome Browser utilities
  * Java environment (>= 1.8)
