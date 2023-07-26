#! /usr/bin/env bash

# Set the number of CPUs/threads for the analysis
export ncpu=1

# Set the Java max memory flag
export java_xmx_flag="-Xmx1g"

# Set the sample ID
export sample_id="truncated_scrnaseq"

# Run Isosceles
singularity exec ../singularity/isosceles.sif Rscript run_isosceles.R

# Run IsoQuant
conda activate isosceles_isoquant
mkdir -p isoquant_results
isoquant.py \
  --threads $ncpu \
  --no_model_construction \
  --data_type nanopore \
  --reference ../reference_data/genome.fasta \
  --genedb ../reference_data/benchmark_transcript_annotations.gtf \
  --bam ../input_data/bam_sim/${sample_id}.bam \
  --read_group tag:BC \
  -o isoquant_results
conda deactivate

# Prepare a FASTQ file from the simulated BAM file for FLAMES
mkdir -p fastq
bedtools bamtofastq -i ../input_data/bam_sim/${sample_id}.bam \
  -fq fastq/${sample_id}.fastq

# Run FLAMES
conda activate isosceles_flames_deps
mkdir -p flames_results
../software/git/FLAMES/python/sc_long_pipeline.py \
  -a ../reference_data/benchmark_transcript_annotations.gtf \
  -i fastq/${sample_id}.fastq \
  -b ../input_data/bam_sim/${sample_id}.bam \
  --outdir flames_results \
  --genomefa ../reference_data/genome.fasta \
  --minimap2_dir $(which minimap2 | xargs dirname)
conda deactivate

# Prepare the reference annotations (refFlat) file for Sicelore
mkdir -p reference_data
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons \
  ../reference_data/benchmark_transcript_annotations.gtf \
  /dev/stdout | \
  awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > \
  reference_data/benchmark_transcript_annotations.refFlat.txt

# Tag the simulated BAM file with gene names for Sicelore
mkdir -p bam_sim
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  AddGeneNameTag \
  I=../input_data/bam_sim/${sample_id}.bam \
  O=bam_sim/${sample_id}_GE.bam \
  REFFLAT=reference_data/benchmark_transcript_annotations.refFlat.txt
samtools index bam_sim/${sample_id}_GE.bam

# Run Sicelore
mkdir -p sicelore_results
java $java_xmx_flag -jar ../software/git/sicelore/Jar/Sicelore-2.0.jar \
  IsoformMatrix \
  I=bam_sim/${sample_id}_GE.bam \
  OUTDIR=sicelore_results \
  PREFIX=sicelore \
  ISOBAM=true \
  GENETAG=GE UMITAG=U8 CELLTAG=BC \
  REFFLAT=reference_data/benchmark_transcript_annotations.refFlat.txt \
  CSV=../reference_data/cell_barcodes.txt \
  DELTA=2 MAXCLIP=150 METHOD=STRICT AMBIGUOUS_ASSIGN=false
