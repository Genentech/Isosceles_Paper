# Mouse E18 brain scRNA-Seq (Lebrigand et al.) data analysis

## Content

  * **download_data.sh** - the script for downloading input data for the analysis
  * **run.sh** - the script for running the analysis scripts
  * **00_run_isosceles.R** - the R script for processing input data using Isosceles (transcript detection and expression quantification)
  * **01_scrnaseq_analysis.R** - the R script for basic scRNA-Seq analysis (dimensionality reduction, clustering etc.)
  * **02_run_dexseq_intra.R** - the R script for running DEXSeq on PSI count data (intra-trajectory analysis)
  * **03_run_dexseq_inter.R** - the R script for running DEXSeq on PSI count data (inter-trajectory analysis)
  * **04_psi_heatmap_intra.R** - the R script for creating the PSI event heatmap (intra-trajectory analysis)
  * **05_psi_heatmap_inter.R** - the R script for creating the PSI event heatmap (inter-trajectory analysis)
  * **06_other_plots.R** - the R script for creating other plots
  * **07_isoswitch_analysis.R** - the R script for isoform switching analysis
  * **results/** - a directory for storing output tabular files and plots

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
