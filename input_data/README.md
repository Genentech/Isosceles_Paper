# Preparing input data for the analysis

## Content

  * **run.sh** - the script for downloading input data
  * **simulate_ont.sh** - a script containing the commands used for simulating Nanopore data
  * **add_bam_tags.py** - a Python script used for adding scRNA-Seq tags to the simulated BAM file (basic benchmarks)
  * **add_bam_tags_ovarian.py** - a Python script used for adding scRNA-Seq tags to simulated BAM files (ovarian cell line analysis)
  * **fastq_ont/** - a directory for storing input Nanopore FASTQ files
  * **fastq_sim/** - a directory for storing simulated Nanopore FASTQ files
  * **bam_sim/** - a directory for storing simulated Nanopore scRNA-Seq BAM files
  * **cellranger/** - a directory for storing Illumina CellRanger output files
  * **read_models/** - a directory for storing selected read model files used for simulating Nanopore data

## Required programs and libraries

  * conda
  * SRA Toolkit
  * samtools (>= 1.10)
