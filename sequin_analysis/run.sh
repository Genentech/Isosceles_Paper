#! /usr/bin/env bash

# Set the number of CPUs/threads for the analysis
export ncpu=1

# Download Sequin annotation data
mkdir -p reference_data
wget -O reference_data/sequin.fasta \
  https://raw.githubusercontent.com/XueyiDong/LongReadRNA/master/sequins/annotations/rnasequin_decoychr_2.4.fa
wget -O reference_data/sequin_annotations.gtf \
  https://raw.githubusercontent.com/XueyiDong/LongReadRNA/master/sequins/annotations/rnasequin_annotation_2.4.gtf
wget -O reference_data/sequin_genes.tsv \
  https://raw.githubusercontent.com/XueyiDong/LongReadRNA/master/sequins/annotations/rnasequin_genes_2.4.tsv
wget -O reference_data/sequin_transcripts.tsv \
  https://raw.githubusercontent.com/XueyiDong/LongReadRNA/master/sequins/annotations/rnasequin_isoforms_2.4.tsv

# Process the Sequin annotations
singularity exec ../singularity/bambu.sif Rscript process_reference_data.R
cat reference_data/sequin_annotations.gtf | awk '$3 == "exon"' > \
  reference_data/sequin_annotations_flames.gtf

# Download Sequin SRA files
mkdir -p sra
wget -O sra/mixA_ont.sra \
  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14286054/SRR14286054
wget -O sra/mixB_ont.sra \
  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14286063/SRR14286063

# Convert Sequin SRA files to FASTQ files
mkdir -p fastq
for sample_id in mixA_ont mixB_ont; do
fasterq-dump -e $ncpu -o fastq/${sample_id}.fastq sra/${sample_id}.sra
gzip fastq/${sample_id}.fastq
done

# Align Nanopore reads to the Sequin sequences with minimap2
mkdir -p bam
for sample_id in mixA_ont mixB_ont; do
../software/bin/minimap2 -t $ncpu \
  -ax splice --secondary=no \
  --junc-bed reference_data/known_introns.bed --junc-bonus 15 \
  reference_data/sequin.fasta \
  fastq/${sample_id}.fastq.gz \
  | samtools sort -o bam/${sample_id}.bam
samtools index bam/${sample_id}.bam
done

# Run Isosceles
singularity exec ../singularity/isosceles.sif Rscript run_isosceles.R

# Run bambu
singularity exec ../singularity/bambu.sif Rscript run_bambu.R

# Run IsoQuant
conda activate isosceles_isoquant
mkdir -p isoquant_results
for sample_id in mixA_ont mixB_ont; do
isoquant.py \
  --threads $ncpu \
  --no_model_construction \
  --data_type nanopore \
  --reference reference_data/sequin.fasta \
  --genedb reference_data/sequin_annotations.gtf \
  --fastq fastq/${sample_id}.fastq.gz \
  -o isoquant_results/${sample_id}
done
conda deactivate

# Run Flair
conda activate isosceles_flair
mkdir -p flair_results
for sample_id in mixA_ont mixB_ont; do
mkdir -p flair_results/${sample_id}
echo sample$'\t'condition$'\t'batch$'\t'fastq/${sample_id}.fastq.gz \
  > flair_results/${sample_id}/read_manifest.tsv
flair 1234 \
  --threads $ncpu \
  --reads fastq/${sample_id}.fastq.gz \
  --genome  reference_data/sequin.fasta \
  --gtf reference_data/sequin_annotations.gtf \
  --junction_bed reference_data/known_introns.bed \
  --reads_manifest flair_results/${sample_id}/read_manifest.tsv \
  --output flair_results/${sample_id}/flair
done
conda deactivate

# Align Nanopore reads to the reference transcriptome with minimap2 for NanoCount
for sample_id in mixA_ont mixB_ont; do
../software/bin/minimap2 -t $ncpu \
  -ax map-ont -p 0 -N 10 \
  reference_data/transcriptome.fasta \
  fastq/${sample_id}.fastq.gz \
  | samtools view -bh > bam/${sample_id}_nanocount.bam
done

# Run NanoCount
conda activate isosceles_nanocount
mkdir -p nanocount_results
for sample_id in mixA_ont mixB_ont; do
NanoCount \
  -i bam/${sample_id}_nanocount.bam \
  -o nanocount_results/${sample_id}.tsv
done
conda deactivate

# Prepare the reference annotations (refgene) files for LIQA
conda activate isosceles_liqa
liqa -task refgene -format gtf \
  -ref reference_data/sequin_annotations.gtf \
  -out reference_data/sequin_annotations.refgene
conda deactivate

# Filter BAM files for LIQA
for sample_id in mixA_ont mixB_ont; do
samtools view \
  bam/${sample_id}.bam \
  -F 2308 -q 50 -O BAM \
  -o bam/${sample_id}_liqa.bam
samtools index bam/${sample_id}_liqa.bam
done

# Run LIQA
conda activate isosceles_liqa
mkdir -p liqa_results
for sample_id in mixA_ont mixB_ont; do
liqa -task quantify \
  -max_distance 10 -f_weight 1 \
  -refgene reference_data/sequin_annotations.refgene \
  -bam bam/${sample_id}_liqa.bam \
  -out liqa_results/${sample_id}.tsv
done
conda deactivate

# Run ESPRESSO_S
conda activate isosceles_espresso_deps
mkdir -p espresso_results
for sample_id in mixA_ont mixB_ont; do
mkdir -p espresso_results/${sample_id}
echo bam/${sample_id}.bam$'\t'sample \
  > espresso_results/${sample_id}/samples.tsv
perl ../software/git/espresso/src/ESPRESSO_S.pl \
  -L espresso_results/${sample_id}/samples.tsv \
  -F reference_data/sequin.fasta \
  -A reference_data/sequin_annotations.gtf \
  -T $ncpu \
  -O espresso_results/${sample_id}
done
conda deactivate

# Run ESPRESSO_C
conda activate isosceles_espresso_deps
for sample_id in mixA_ont mixB_ont; do
perl ../software/git/espresso/src/ESPRESSO_C.pl \
  -I espresso_results/${sample_id} \
  -F reference_data/sequin.fasta \
  -T $ncpu \
  -X 0
done
conda deactivate

# Run ESPRESSO_Q
conda activate isosceles_espresso_deps
for sample_id in mixA_ont mixB_ont; do
perl ../software/git/espresso/src/ESPRESSO_Q.pl \
  -L espresso_results/${sample_id}/samples.tsv.updated \
  -A reference_data/sequin_annotations.gtf \
  -T $ncpu
done
conda deactivate
