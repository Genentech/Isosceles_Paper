import itertools
import pysam

cell_barcode = "AAACCCAAGTATTGCC"

for sample_id in ["SK-OV-3",  "IGROV-1", "OVMANA", "OVKATE", "OVTOKO", "COV362"]:
    umi_iterator = itertools.product("ACGT", repeat=12)
    samfile_path = "bam/{0}_Rep2.bam".format(sample_id)
    outfile_path = "bam/{0}_Rep2_sicelore.bam".format(sample_id)
    with pysam.AlignmentFile(samfile_path, "rb") as samfile:
        with pysam.AlignmentFile(outfile_path, "wb", template=samfile) as outfile:
            for read in samfile:
                umi_sequence = "".join(next(umi_iterator))
                read.set_tag("BC", cell_barcode)
                read.set_tag("U8", umi_sequence)
                outfile.write(read)
