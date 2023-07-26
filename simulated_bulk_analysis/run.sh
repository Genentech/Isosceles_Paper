#! /usr/bin/env bash

# Set the number of CPUs/threads for the analysis
export ncpu=2

# Set the sample ID
export sample_id="truncated_bulk_rnaseq"

# Align Nanopore reads to the reference genome with minimap2
mkdir -p bam
../software/bin/minimap2 -t $ncpu \
  -ax splice --secondary=no \
  --junc-bed ../reference_data/known_introns.bed --junc-bonus 15 \
  ../reference_data/genome.fasta \
  ../input_data/fastq_sim/${sample_id}.fastq.gz \
  | samtools sort -o bam/${sample_id}.bam
samtools index bam/${sample_id}.bam

# Run StringTie
conda activate isosceles_stringtie
mkdir -p stringtie_results
for perc_down in 10 20 30; do
stringtie -L -p $ncpu \
  -G ../reference_data/benchmark_downsampled_${perc_down}.gtf \
  bam/${sample_id}.bam \
  | awk '$7 != "."' \
  > stringtie_results/${sample_id}_denovo_${perc_down}_raw.gtf
stringtie --merge \
  -G ../reference_data/benchmark_downsampled_${perc_down}.gtf \
  -o stringtie_results/${sample_id}_denovo_${perc_down}.gtf \
  stringtie_results/${sample_id}_denovo_${perc_down}_raw.gtf
cat stringtie_results/${sample_id}_denovo_${perc_down}.gtf | \
  awk '$3 == "exon"' > \
  stringtie_results/${sample_id}_denovo_${perc_down}_exon.gtf
done
conda deactivate

# Extract StringTie transcript sequences
singularity exec ../singularity/bambu.sif Rscript prepare_stringtie_fasta.R

# Run IsoQuant (transcript quantification)
conda activate isosceles_isoquant
mkdir -p isoquant_results
isoquant.py \
  --threads $ncpu \
  --no_model_construction \
  --data_type nanopore \
  --reference ../reference_data/genome.fasta \
  --genedb ../reference_data/benchmark_transcript_annotations.gtf \
  --bam bam/${sample_id}.bam \
  -o isoquant_results/${sample_id}_quant
conda deactivate

# Run IsoQuant (de novo detection)
conda activate isosceles_isoquant
mkdir -p isoquant_results
for perc_down in 10 20 30; do
isoquant.py \
  --threads $ncpu \
  --data_type nanopore \
  --reference ../reference_data/genome.fasta \
  --genedb ../reference_data/benchmark_downsampled_${perc_down}.gtf \
  --bam bam/${sample_id}.bam \
  -o isoquant_results/${sample_id}_denovo_${perc_down}
cat ../reference_data/benchmark_downsampled_${perc_down}.gtf > \
  isoquant_results/${sample_id}_denovo_${perc_down}/extended_annotations.gtf
cat isoquant_results/${sample_id}_denovo_${perc_down}/00_${sample_id}/00_${sample_id}.transcript_models.gtf \
  | grep "nic" >> \
  isoquant_results/${sample_id}_denovo_${perc_down}/extended_annotations.gtf
isoquant.py \
  --threads $ncpu \
  --no_model_construction \
  --data_type nanopore \
  --reference ../reference_data/genome.fasta \
  --genedb isoquant_results/${sample_id}_denovo_${perc_down}/extended_annotations.gtf \
  --bam bam/${sample_id}.bam \
  -o isoquant_results/${sample_id}_quant_denovo_${perc_down}
done
conda deactivate

# Run IsoQuant (de novo detection using StringTie)
conda activate isosceles_isoquant
mkdir -p isoquant_results
for perc_down in 10 20 30; do
isoquant.py \
  --threads $ncpu \
  --data_type nanopore \
  --reference ../reference_data/genome.fasta \
  --genedb stringtie_results/${sample_id}_denovo_${perc_down}.gtf \
  --bam bam/${sample_id}.bam \
  -o isoquant_results/${sample_id}_stringtie_${perc_down}
cat stringtie_results/${sample_id}_denovo_${perc_down}.gtf > \
  isoquant_results/${sample_id}_stringtie_${perc_down}/extended_annotations.gtf
cat isoquant_results/${sample_id}_stringtie_${perc_down}/00_${sample_id}/00_${sample_id}.transcript_models.gtf \
  | grep "nic" >> \
  isoquant_results/${sample_id}_stringtie_${perc_down}/extended_annotations.gtf
isoquant.py \
  --threads $ncpu \
  --no_model_construction \
  --data_type nanopore \
  --reference ../reference_data/genome.fasta \
  --genedb isoquant_results/${sample_id}_stringtie_${perc_down}/extended_annotations.gtf \
  --bam bam/${sample_id}.bam \
  -o isoquant_results/${sample_id}_quant_stringtie_${perc_down}
done
conda deactivate

# Run Isosceles
singularity exec ../singularity/isosceles.sif Rscript run_isosceles.R

# Run bambu
singularity exec ../singularity/bambu.sif Rscript run_bambu.R

# Run Flair (transcript quantification)
conda activate isosceles_flair
mkdir -p flair_results/${sample_id}_quant
echo sample$'\t'condition$'\t'batch$'\t'../input_data/fastq_sim/${sample_id}.fastq.gz \
  > flair_results/${sample_id}_quant/read_manifest.tsv
flair 1234 \
  --threads $ncpu \
  --reads ../input_data/fastq_sim/${sample_id}.fastq.gz \
  --genome  ../reference_data/genome.fasta \
  --gtf ../reference_data/benchmark_transcript_annotations.gtf \
  --junction_bed ../reference_data/known_introns.bed \
  --reads_manifest flair_results/${sample_id}_quant/read_manifest.tsv \
  --output flair_results/${sample_id}_quant/flair
conda deactivate

# Run Flair (de novo detection)
conda activate isosceles_flair
for perc_down in 10 20 30; do
mkdir -p flair_results/${sample_id}_denovo_${perc_down}
echo sample$'\t'condition$'\t'batch$'\t'../input_data/fastq_sim/${sample_id}.fastq.gz \
  > flair_results/${sample_id}_denovo_${perc_down}/read_manifest.tsv
flair 1234 \
  --threads $ncpu \
  --reads ../input_data/fastq_sim/${sample_id}.fastq.gz \
  --genome  ../reference_data/genome.fasta \
  --gtf ../reference_data/benchmark_downsampled_${perc_down}.gtf \
  --junction_bed ../reference_data/known_introns.bed \
  --reads_manifest flair_results/${sample_id}_denovo_${perc_down}/read_manifest.tsv \
  --output flair_results/${sample_id}_denovo_${perc_down}/flair
done
conda deactivate

# Run Flair (de novo detection using StringTie)
conda activate isosceles_flair
for perc_down in 10 20 30; do
mkdir -p flair_results/${sample_id}_stringtie_${perc_down}
echo sample$'\t'condition$'\t'batch$'\t'../input_data/fastq_sim/${sample_id}.fastq.gz \
  > flair_results/${sample_id}_stringtie_${perc_down}/read_manifest.tsv
flair 1234 \
  --threads $ncpu \
  --reads ../input_data/fastq_sim/${sample_id}.fastq.gz \
  --genome  ../reference_data/genome.fasta \
  --gtf stringtie_results/${sample_id}_denovo_${perc_down}.gtf \
  --junction_bed ../reference_data/known_introns.bed \
  --reads_manifest flair_results/${sample_id}_stringtie_${perc_down}/read_manifest.tsv \
  --output flair_results/${sample_id}_stringtie_${perc_down}/flair
done
conda deactivate

# Align Nanopore reads to the reference transcriptome with minimap2 for NanoCount
../software/bin/minimap2 -t $ncpu \
  -ax map-ont -p 0 -N 10 \
  ../reference_data/benchmark_transcript_sequences.fasta \
  ../input_data/fastq_sim/${sample_id}.fastq.gz \
  | samtools view -bh > bam/${sample_id}_nanocount_quant.bam

# Run NanoCount (transcript quantification)
conda activate isosceles_nanocount
mkdir -p nanocount_results
NanoCount \
  -i bam/${sample_id}_nanocount_quant.bam \
  -o nanocount_results/${sample_id}_quant.tsv
conda deactivate

# Align Nanopore reads to StringTie transcriptomes with minimap2 for NanoCount
for perc_down in 10 20 30; do
../software/bin/minimap2 -t $ncpu \
  -ax map-ont -p 0 -N 10 \
  stringtie_results/${sample_id}_denovo_${perc_down}.fasta \
  ../input_data/fastq_sim/${sample_id}.fastq.gz \
  | samtools view -bh > bam/${sample_id}_nanocount_stringtie_${perc_down}.bam
done

# Run NanoCount (de novo detection using StringTie)
conda activate isosceles_nanocount
mkdir -p nanocount_results
for perc_down in 10 20 30; do
NanoCount \
  -i bam/${sample_id}_nanocount_stringtie_${perc_down}.bam \
  -o nanocount_results/${sample_id}_stringtie_${perc_down}.tsv
done
conda deactivate

# Prepare the reference annotations (refgene) files for LIQA
conda activate isosceles_liqa
mkdir -p reference_data
liqa -task refgene -format gtf \
  -ref ../reference_data/benchmark_transcript_annotations.gtf \
  -out reference_data/benchmark_transcript_annotations.refgene
for perc_down in 10 20 30; do
liqa -task refgene -format gtf \
  -ref stringtie_results/${sample_id}_denovo_${perc_down}.gtf \
  -out stringtie_results/${sample_id}_denovo_${perc_down}.refgene
done
conda deactivate

# Filter the BAM file for LIQA
samtools view \
  bam/${sample_id}.bam \
  -F 2308 -q 50 -O BAM \
  -o bam/${sample_id}_liqa.bam
samtools index bam/${sample_id}_liqa.bam

# Run LIQA (transcript quantification)
conda activate isosceles_liqa
mkdir -p liqa_results
liqa -task quantify \
  -max_distance 10 -f_weight 1 \
  -refgene reference_data/benchmark_transcript_annotations.refgene \
  -bam bam/${sample_id}_liqa.bam \
  -out liqa_results/${sample_id}_quant.tsv
conda deactivate

# Run LIQA (de novo detection using StringTie)
conda activate isosceles_liqa
mkdir -p liqa_results
for perc_down in 10 20 30; do
liqa -task quantify \
  -max_distance 10 -f_weight 1 \
  -refgene stringtie_results/${sample_id}_denovo_${perc_down}.refgene \
  -bam bam/${sample_id}_liqa.bam \
  -out liqa_results/${sample_id}_stringtie_${perc_down}.tsv
done
conda deactivate

# Run ESPRESSO_S (transcript quantification)
conda activate isosceles_espresso_deps
mkdir -p espresso_results/${sample_id}_quant
echo bam/${sample_id}.bam$'\t'sample \
  > espresso_results/${sample_id}_quant/samples.tsv
perl ../software/git/espresso/src/ESPRESSO_S.pl \
  -L espresso_results/${sample_id}_quant/samples.tsv \
  -F ../reference_data/genome.fasta \
  -A ../reference_data/benchmark_transcript_annotations.gtf \
  -M MT \
  -T $ncpu \
  -O espresso_results/${sample_id}_quant
conda deactivate

# Run ESPRESSO_C (transcript quantification)
conda activate isosceles_espresso_deps
perl ../software/git/espresso/src/ESPRESSO_C.pl \
  -I espresso_results/${sample_id}_quant \
  -F ../reference_data/genome.fasta \
  -T $ncpu \
  -X 0
conda deactivate

# Run ESPRESSO_Q (transcript quantification)
conda activate isosceles_espresso_deps
perl ../software/git/espresso/src/ESPRESSO_Q.pl \
  -L espresso_results/${sample_id}_quant/samples.tsv.updated \
  -A ../reference_data/benchmark_transcript_annotations.gtf \
  -T $ncpu
conda deactivate

# Run ESPRESSO_S (de novo detection)
conda activate isosceles_espresso_deps
for perc_down in 10 20 30; do
mkdir -p espresso_results/${sample_id}_denovo_${perc_down}
echo bam/${sample_id}.bam$'\t'sample \
  > espresso_results/${sample_id}_denovo_${perc_down}/samples.tsv
perl ../software/git/espresso/src/ESPRESSO_S.pl \
  -L espresso_results/${sample_id}_denovo_${perc_down}/samples.tsv \
  -F ../reference_data/genome.fasta \
  -A ../reference_data/benchmark_downsampled_${perc_down}.gtf \
  -M MT \
  -T $ncpu \
  -O espresso_results/${sample_id}_denovo_${perc_down}
done
conda deactivate

# Run ESPRESSO_C (de novo detection)
conda activate isosceles_espresso_deps
for perc_down in 10 20 30; do
perl ../software/git/espresso/src/ESPRESSO_C.pl \
  -I espresso_results/${sample_id}_denovo_${perc_down} \
  -F ../reference_data/genome.fasta \
  -T $ncpu \
  -X 0
done
conda deactivate

# Run ESPRESSO_Q (de novo detection)
conda activate isosceles_espresso_deps
for perc_down in 10 20 30; do
perl ../software/git/espresso/src/ESPRESSO_Q.pl \
  -L espresso_results/${sample_id}_denovo_${perc_down}/samples.tsv.updated \
  -A ../reference_data/benchmark_downsampled_${perc_down}.gtf \
  -T $ncpu
done
conda deactivate

# Run ESPRESSO_S (de novo detection using StringTie)
conda activate isosceles_espresso_deps
for perc_down in 10 20 30; do
mkdir -p espresso_results/${sample_id}_stringtie_${perc_down}
echo bam/${sample_id}.bam$'\t'sample \
  > espresso_results/${sample_id}_stringtie_${perc_down}/samples.tsv
perl ../software/git/espresso/src/ESPRESSO_S.pl \
  -L espresso_results/${sample_id}_stringtie_${perc_down}/samples.tsv \
  -F ../reference_data/genome.fasta \
  -A stringtie_results/${sample_id}_denovo_${perc_down}_exon.gtf \
  -M MT \
  -T $ncpu \
  -O espresso_results/${sample_id}_stringtie_${perc_down}
done
conda deactivate

# Run ESPRESSO_C (de novo detection using StringTie)
conda activate isosceles_espresso_deps
for perc_down in 10 20 30; do
perl ../software/git/espresso/src/ESPRESSO_C.pl \
  -I espresso_results/${sample_id}_stringtie_${perc_down} \
  -F ../reference_data/genome.fasta \
  -T $ncpu \
  -X 0
done
conda deactivate

# Run ESPRESSO_Q (de novo detection using StringTie)
conda activate isosceles_espresso_deps
for perc_down in 10 20 30; do
perl ../software/git/espresso/src/ESPRESSO_Q.pl \
  -L espresso_results/${sample_id}_stringtie_${perc_down}/samples.tsv.updated \
  -A stringtie_results/${sample_id}_denovo_${perc_down}_exon.gtf \
  -T $ncpu
done
conda deactivate
