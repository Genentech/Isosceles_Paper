import pysam

cell_barcodes = []
with open("../reference_data/cell_barcodes.txt") as infile:
    for line in infile:
        cell_barcodes.append(line.strip())

umi_sequences = []
with open("../reference_data/umi_sequences.txt") as infile:
    for line in infile:
        umi_sequences.append(line.strip())

for i in range(100):
    cell_barcode = cell_barcodes[i]
    samfile_path = "bam/subsample/truncated_scrnaseq_{0}.bam".format(i + 1)
    outfile_path = "bam/subsample_tags/truncated_scrnaseq_{0}.bam".format(i + 1)
    with pysam.AlignmentFile(samfile_path, "rb") as samfile:
        with pysam.AlignmentFile(outfile_path, "wb", template=samfile) as outfile:
            for j, read in enumerate(samfile):
                umi_sequence = umi_sequences[j]
                read.set_tag("BC", cell_barcode)
                read.set_tag("U8", umi_sequence)
                read.query_name = "{0}_{1}#{2}".format(cell_barcode, umi_sequence, read.query_name)
                outfile.write(read)
