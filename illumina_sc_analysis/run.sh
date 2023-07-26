#! /usr/bin/env bash

# Set the number of CPUs/threads for the analysis
export ncpu=1

# Set the Java max memory flag
export java_xmx_flag="-Xmx1g"

# Run souporcell
singularity exec \
  shub://wheaton5/souporcell \
  souporcell_pipeline.py \
  -i ../input_data/cellranger/possorted_genome_bam.bam \
  -b ../input_data/cellranger/barcodes.tsv \
  -f ../reference_data/genome.fasta \
  -t $ncpu \
  -k 3 \
  -o souporcell
cp -f souporcell/clusters.tsv data/souporcell_clusters.tsv

# Parse Illumina data using Sicelore
java $java_xmx_flag -jar ../software/git/sicelore/Jar/IlluminaParser-1.0.jar \
  --inFileIllumina ../input_data/cellranger/possorted_genome_bam.bam \
  --tsv ../input_data/cellranger/barcodes.tsv \
  --cellBCflag CB --umiFlag UB --geneFlag GX \
  --outFile data/illumina_parsed.obj

# Make Illumina barcodes usable with Nanopore data
cat ../input_data/cellranger/barcodes.tsv | sed 's/-1//g' > data/barcodes.csv
