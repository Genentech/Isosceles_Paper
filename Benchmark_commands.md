# Benchmark command list

## NanoSim 

### Building read models
```bash
read_analysis.py transcriptome -t {threads} --no_intron_retention -i {fastq_file} -rt {transcriptome_fasta_file} -rg  {genome_fasta_file} -o {read_model_file_prefix}
```

### Simulating Nanopore reads
```bash
simulator.py transcriptome -t {threads} --fastq --no_model_ir -b albacore -r cDNA_1D -n 100000000 -c {read_model_file_prefix} -rt {transcriptome_fasta_file} -e {expression_values_tabular_file} -o {simulated_reads_file_prefix}
cat {simulated_reads_file_prefix}_aligned_reads.fastq | head -n 48000000 | gzip -c > {fastq_file}
```

## minimap2

### Aligning Nanopore reads
```bash
minimap2 -t {threads} -ax splice --secondary=no --junc-bed {junction_bed_file} --junc-bonus 15 {genome_fasta_file} {fastq_file} | samtools sort -o {bam_file}
samtools index {bam_file}
```

### Aligning Pacbio reads
```bash
minimap2 -t {threads} -ax splice:hq --secondary=no --junc-bed {junction_bed_file} --junc-bonus 15 {genome_fasta_file} {fastq_file} | samtools sort -o {bam_file}
samtools index {bam_file}
```

## Simulating scRNA-Seq data

### Subsampling the aligned simulated reads
```bash
samtools view -b -s {subsampling_seed}.{subsampling_fraction} -o {subsampled_bam_file} {input_bam_file}
```

### <a id="add_tags" /> Adding cell barcode and UMI tags
```python
bc_iterator = itertools.product("ACGT", repeat=16)
umi_iterator = itertools.product("ACGT", repeat=12)
with pysam.AlignmentFile({subsampled_bam_file}, "rb") as samfile:
    with pysam.AlignmentFile({tagged_subsampled_bam_file}, "wb", template=samfile) as outfile:
        cell_barcode = "".join(next(bc_iterator))
        for read in samfile:
            umi_sequence = "".join(next(umi_iterator))
            read.set_tag("BC", cell_barcode)
            read.set_tag("U8", umi_sequence)
            read.query_name = "{0}_{1}#{2}".format(cell_barcode, umi_sequence, read.query_name)
            outfile.write(read)
```

### Merging and sorting the subsampled BAM files
```bash
samtools merge {unsorted_bam_file} {tagged_subsampled_bam_file_directory}/*.bam
samtools sort -o {bam_file} {unsorted_bam_file}
samtools index {bam_file}
```

### Preparing a FASTQ file from the simulated BAM file
```bash
bedtools bamtofastq -i {bam_file} -fq {fastq_file}
```

## StringTie

### Running StringTie
```bash
stringtie -L -p {threads} -G {gtf_file} {bam_file} | awk '\$7 != "."' > {stringtie_raw_gtf_file}
stringtie --merge -G {gtf_file} -o {stringtie_gtf_file} {stringtie_raw_gtf_file}
```

### Extract StringTie transcriptome sequences
```r
genome_seq <- Biostrings::readDNAStringSet({genome_fasta}, format = "fasta")
names(genome_seq) <- sapply(strsplit(names(genome_seq), "\\s+"), "[", 1)
txdb <- GenomicFeatures::makeTxDbFromGFF({stringtie_gtf_file})
transcriptome_seq <- GenomicFeatures::extractTranscriptSeqs(genome_seq, txdb, use.names = TRUE)
Biostrings::writeXStringSet(transcriptome_seq, {stringtie_transcriptome_fasta_file}, format = "fasta")
```

## bambu 

### Transcript quantification
```r
txdb <- GenomicFeatures::makeTxDbFromGFF({gtf_file})
annotations <- bambu::prepareAnnotations(txdb)
se <- bambu::bambu(reads = {bam_file}, annotations = annotations, genome = {genome_fasta_file}, discovery = FALSE, ncore = {threads}, yieldSize = 1e6, lowMemory = TRUE)
bambu::writeBambuOutput(se, path = {output_directory}, prefix = {output_file_prefix})
```

### De novo detection
```r
txdb <- GenomicFeatures::makeTxDbFromGFF({gtf_file})
annotations <- bambu::prepareAnnotations(txdb)
se <- bambu::bambu(reads = {bam_file}, annotations = annotations, genome = {genome_fasta_file}, discovery = TRUE, ncore = {threads}, yieldSize = 1e6, lowMemory = TRUE)
bambu::writeBambuOutput(se, path = {output_directory}, prefix = {output_file_prefix})
```

## Flair
```bash
echo sample$'\t'condition$'\t'batch$'\t'{fastq_file} > {read_manifest_file}
flair 1234 --threads {threads} --reads {fastq_file} --genome {genome_fasta_file} --gtf {gtf_file} --junction_bed {junction_bed_file --reads_manifest {read_manifest_file} --output {output_file_prefix}
```

## Nanocount
```bash
minimap2 -t {threads} -ax map-ont -p 0 -N 10 {transcriptome_fasta_file} {fastq_file} | samtools view -bh > {bam_file}
NanoCount -i {bam_file} -o {output_file_path}
```

## LIQA
```bash
liqa -task refgene -format gtf -ref {gtf_file} -out {refgene_file}
samtools view {bam_file} -F 2308 -q 50 -O BAM -o {filtered_bam_file}
liqa -task quantify -max_distance 10 -f_weight 1 -refgene {refgene_file} -bam {filtered_bam_file} -out {output_file_path}
```

## ESPRESSO
```bash
echo {bam_file}$'\t'sample > {sample_table_file}
perl ESPRESSO_S.pl -L {sample_table_file} -F {genome_fasta_file} -A {gtf_file} -M MT -T {threads} -O {output_directory}
perl ESPRESSO_C.pl -I {output_directory} -F {genome_fasta_file} -T {threads} -X 0
perl ESPRESSO_Q.pl -L {output_directory}/samples.tsv.updated -A {gtf_file} -T {threads}
```

## souporcell 
```bash
singularity exec shub://wheaton5/souporcell souporcell_pipeline.py -i {cellranger_directory}/possorted_genome_bam.bam -b {cell_barcodes_file} -f {genome_fasta_file} -t {threads} -k 3 -o {output_directory}
```

## Sicelore 

### <a id="sicelore_preprocessing" /> Data preprocessing
  * Sicelore workflow is described here: https://github.com/ucagenomix/sicelore
  * For non-simulated scRNA-Seq data, steps 1-8 need to be run
    * For step 7, we used the flags shown in the [minimap2](#minimap2) section
  * For simulated scRNA-Seq data, only step 8 (AddGeneNameTag) needs to be run
  * For bulk RNA-Seq data we need two steps:
    * aligned reads need to be processed using the [Step 2](#add_tags) code of the 'Simulating scRNA-Seq data' section
    * step 8 (AddGeneNameTag) of the Sicelore workflow needs to be run

### Preparing the refFlat file
```bash
gtfToGenePred -genePredExt -ignoreGroupsWithoutExons {gtf_file} /dev/stdout | awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' > {refflat_file}
```

### Running Sicelore
```bash
java -jar Sicelore-2.0.jar IsoformMatrix I={bam_file} OUTDIR={output_directory} PREFIX={output_file_prefix} ISOBAM=true GENETAG=GE UMITAG=U8 CELLTAG=BC REFFLAT={refflat_file} CSV={cell_barcodes_file} DELTA=2 MAXCLIP=150 METHOD=STRICT AMBIGUOUS_ASSIGN=false
```

## wf-single-cell 

### Preparing Cell Ranger custom reference data
```bash
cellranger mkref --nthreads={threads} --genome={cellranger_directory} --fasta={genome_fasta_file} --genes={gtf_file}
gunzip {cellranger_directory}/genes/genes.gtf.gz
```

### Running the wf-single-cell workflow
```bash
export TMPDIR={temp_directory}
nextflow run epi2me-labs/wf-single-cell -r v1.1.0 -process.executor='local' -c {config_file} -profile singularity -work-dir {work_directory} --max_threads {threads} --resources_mm2_max_threads {threads} --resources_mm2_flags="--junc-bonus 15" --merge_bam True --fastq {fastq_file} --kit_name 3prime --kit_version v3 --expected_cells 2000 --ref_genome_dir {cellranger_directory} --plot_umaps --out_dir {output_directory}
```

### Filtering alignments with the correct UMI length
```bash
samtools view -b --expr 'length([UB]) == 12' {input_bam_file} > {bam_file}
samtools index {bam_file}
```

### Deduplicating reads using UMI and mapping coordinates
```bash
umi_tools dedup --extract-umi-method=tag --method=directional --per-gene --per-cell --umi-tag=UB --cell-tag=CB --gene-tag=GN --stdin {input_bam_file} --output-stats={stats_file_prefix} --stdout {bam_file}
samtools index {bam_file}
```

## IsoQuant

### Simulated bulk RNA-Seq data
  * For simulated bulk RNA-Seq data BAM files (rather than FASTQ files) are used as input, to account for the fact that all other programs have access to alignments using the whole known intron set

#### Transcript quantification
```bash
isoquant.py --threads {threads} --no_model_construction --data_type nanopore --reference {genome_fasta_file} --genedb {gtf_file} --bam {bam_file} -o {output_directory}
```

#### De novo detection
  * We need to run IsoQuant twice (once in de novo detection detection mode, and then just for transcript quantification with extended annotations), as IsoQuant doesn't seem to provide transcript expression values for all the known and discovered transcripts by default
```bash
isoquant.py --threads {threads} --data_type nanopore --reference {genome_fasta_file} --genedb {gtf_file} --bam {bam_file} -o {output_directory}
cat {gtf_file} > {output_directory}/extended_annotations.gtf
cat {output_directory}/00_{bam_file_basename}/00_{bam_file_basename}.transcript_models.gtf | grep "nic" >> {output_directory}/extended_annotations.gtf
isoquant.py --threads {threads} --no_model_construction --data_type nanopore --reference {genome_fasta_file} --genedb {output_directory}/extended_annotations.gtf --bam {bam_file} -o {output_directory_2}
```

### Non-simulated bulk RNA-Seq data
```bash
isoquant.py --threads {threads} --no_model_construction --data_type nanopore --reference {genome_fasta_file} --genedb {gtf_file} --fastq {fastq_file} -o {output_directory}
```

### scRNA-Seq data
  * For non-simulated scRNA-Seq data, IsoQuant requires a BAM file prepared using the [Sicelore preprocessing workflow](#sicelore_preprocessing) as input
```bash
isoquant.py --threads {threads} --no_model_construction --data_type nanopore --reference {genome_fasta_file} --genedb {gtf_file} --bam {bam_file} --read_group tag:BC -o {output_directory}
```

## FLAMES 

### Bulk RNA-Seq data
  * For bulk RNA-Seq data, we used the default config file (config_sclr_nanopore_bulk.json) as a template, changing the 'has_UMI' value to false, and the 'strand_specific' value to 0
```bash
bulk_long_pipeline.py --config_file {config_file} -a {gtf_file} -i {fastq_file} --outdir {output_directory} --genomefa {genome_fasta_file} --minimap2_dir {minimap2_dir}
```

### Simulated scRNA-Seq data
```bash
sc_long_pipeline.py -a {gtf_file} -i {fastq_file} -b {bam_file} --outdir {output_directory} --genomefa {genome_fasta_file} --minimap2_dir {minimap2_dir}
```

### Non-simulated scRNA-Seq data
```bash
match_cell_barcode {input_fastq_directory} {barcode_stats_file} {fastq_file} {cell_barcodes_file} 1 12
sc_long_pipeline.py -a {gtf_file} -i {fastq_file} --outdir {output_directory} --genomefa {genome_fasta_file} --minimap2_dir {minimap2_dir}
```
