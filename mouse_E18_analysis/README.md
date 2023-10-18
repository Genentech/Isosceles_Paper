# Mouse E18 scRNA-Seq (Lebrigand et al.) data analysis

## Content

  * **download_data.sh** - the script for downloading input data for the analysis
  * **00_run_isosceles.R** - the R script for processing input data using Isosceles (transcript detection and expression quantification)

## Input files to download

  * **genome.fasta** - a FASTA file containing reference genome sequences (GRCm38)
  * **gencode_m25.gtf** - a GTF file containing reference genome annotations (source: GENCODE M25)
  * **known_introns.bed** - a BED file containing known intron positions (GENCODE M25 and VastDB mm10 introns)
  * **190_cells.bam** - a BAM file created using the Sicelore workflow (190 cells sequencing)
  * **190_cells.bam.bai** - an index for the 190_cells.bam file
  * **951_cells.bam** - a BAM file created using the Sicelore workflow (951 cells sequencing)
  * **951_cells.bam.bai** - an index for the 951_cells.bam file

## Required programs and libraries

  * Singularity (>= 3.5.1)
