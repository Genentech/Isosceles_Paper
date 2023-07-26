# Preparing reference genome and annotations data for the analysis

## Content

  * **run.sh** - the script for downloading the reference data
  
## Files to download

  * **genome.fasta** - a FASTA file containing reference genome sequences (GRCh38)
  * **transcriptome.fasta** - a FASTA file containing reference transcriptome sequences (source: Ensembl 90)
  * **ensembl_90.gtf** - a GTF file containing reference genome annotations (source: Ensembl 90)
  * **known_introns.bed** - a BED file containing annotated intron positions
  * **cell_barcodes.txt** - a text file containing unique cell barcode sequences
  * **umi_sequences.txt** - a text file containing unique UMI sequences
  * **benchmark_transcript_sequences.fasta** - a FASTA file containing simulated transcript sequences
  * **benchmark_transcript_annotations.gtf** - a GTF file containing simulated transcript annotations
  * **benchmark_transcript_expression.tab** - a tabular file containing simulated transcript expression values
  * **benchmark_transcript_gene.tab** - a tabular file containing simulated transcript to gene mapping data
  * **benchmark_downsampled_\*.gtf** - GTF files containing downsampled (10%, 20% and 30%) simulated transcript annotations
  * **benchmark_downsampled_\*\_kept.tab** - text files containing kept transcript IDs for downsampled (10%, 20% and 30%) simulated transcript data
  * **benchmark_downsampled_\*\_dropped.tab** - text files containing dropped transcript IDs for downsampled (10%, 20% and 30%) simulated transcript data
