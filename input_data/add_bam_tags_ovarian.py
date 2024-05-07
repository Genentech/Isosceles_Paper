import itertools
import random
import pysam

cell_barcodes = []
with open("../reference_data/cell_barcodes.txt") as infile:
    for line in infile:
        cell_barcodes.append(line.strip())

for sample_id in ["SK-OV-3",  "IGROV-1", "OVMANA", "OVKATE", "OVTOKO", "COV362"]:
    umi_iterator = itertools.product("ACGT", repeat=12)
    samfile_path = "bam/{0}_Rep1.bam".format(sample_id)
    outfile_path = "bam_sim/{0}_Rep1.bam".format(sample_id)
    with pysam.AlignmentFile(samfile_path, "rb") as samfile:
        with pysam.AlignmentFile(outfile_path, "wb", template=samfile) as outfile:
            for read in samfile:
                cell_barcode = random.choice(cell_barcodes)
                umi_sequence = "".join(next(umi_iterator))
                read.set_tag("BC", cell_barcode)
                read.set_tag("U8", umi_sequence)
                read.query_name = "{0}_{1}#{2}".format(cell_barcode, umi_sequence, read.query_name)
                outfile.write(read)
