#! /usr/bin/env bash

# Set the number of CPUs/threads for the analysis
export ncpu=1

# Build the read model from bulk RNA-Seq ONT data (mix of three cell lines)
conda activate isosceles_nanosim
mkdir -p read_models/bulk_rnaseq
read_analysis.py transcriptome \
  -t $ncpu --no_intron_retention \
  -i <( zcat fastq_ont/LIB5432312_SAM24385455.fastq.gz fastq_ont/LIB5432314_SAM24385457.fastq.gz fastq_ont/LIB5432315_SAM24385458.fastq.gz ) \
  -rt ../reference_data/transcriptome.fasta \
  -rg  ../reference_data/genome.fasta \
  -o read_models/bulk_rnaseq/training
conda deactivate

# Build the read model from scRNA-Seq ONT data
conda activate isosceles_nanosim
mkdir -p read_models/scrnaseq
read_analysis.py transcriptome \
  -t $ncpu --no_intron_retention \
  -i <( zcat fastq_ont/LIB5445493_SAM24404003.fastq.gz ) \
  -rt ../reference_data/transcriptome.fasta \
  -rg  ../reference_data/genome.fasta \
  -o read_models/scrnaseq/training
conda deactivate

# Simulate reads (bulk RNA-Seq error model)
conda activate isosceles_nanosim
mkdir -p nanosim_results
simulator.py transcriptome \
  -t $ncpu --fastq --no_model_ir -b albacore -r cDNA_1D \
  -n 100000000 \
  -c read_models/bulk_rnaseq/training \
  -rt ../reference_data/benchmark_transcript_sequences.fasta \
  -e ../reference_data/benchmark_transcript_expression.tab \
  -o nanosim_results/truncated_bulk_rnaseq
conda deactivate

# Simulate reads (scRNA-Seq error model)
conda activate isosceles_nanosim
mkdir -p nanosim_results
simulator.py transcriptome \
  -t $ncpu --fastq --no_model_ir -b albacore -r cDNA_1D \
  -n 100000000 \
  -c read_models/scrnaseq/training \
  -rt ../reference_data/benchmark_transcript_sequences.fasta \
  -e ../reference_data/benchmark_transcript_expression.tab \
  -o nanosim_results/truncated_scrnaseq
conda deactivate

# Select the first 12 million aligned reads from the simulated data
cat nanosim_results/truncated_bulk_rnaseq_aligned_reads.fastq | \
  head -n 48000000 | gzip -c > \
  fastq_sim/truncated_bulk_rnaseq.fastq.gz
cat nanosim_results/truncated_scrnaseq_aligned_reads.fastq | \
  head -n 48000000 | gzip -c > \
  fastq_sim/truncated_scrnaseq.fastq.gz

# Simulate the Nanopore scRNA-Seq BAM file
## Align Nanopore reads to the reference genome with minimap2
mkdir -p bam
../software/bin/minimap2 -t $ncpu \
  -ax splice --secondary=no \
  --junc-bed ../reference_data/known_introns.bed --junc-bonus 15 \
  ../reference_data/genome.fasta \
  fastq_sim/truncated_scrnaseq.fastq.gz \
  | samtools sort -o bam/truncated_scrnaseq.bam
samtools index bam/truncated_scrnaseq.bam
## Subsample the BAM file containing aligned simulated reads
mkdir -p bam/subsample
for i in $(seq 1 100); do
samtools view -b -s ${i}.000833 \
  -o bam/subsample/truncated_scrnaseq_${i}.bam \
  bam/truncated_scrnaseq.bam
done
## Add cell barcode and UMI tags to the subsampled BAM files
mkdir -p bam/subsample_tags
conda activate isosceles_nanocount
python add_bam_tags.py
conda deactivate
## Merge and sort the subsampled BAM files
samtools merge bam_sim/truncated_scrnaseq_unsorted.bam \
  bam/subsample_tags/*.bam
samtools sort -o bam_sim/truncated_scrnaseq.bam \
  bam_sim/truncated_scrnaseq_unsorted.bam
samtools index bam_sim/truncated_scrnaseq.bam
rm -f bam_sim/truncated_scrnaseq_unsorted.bam
