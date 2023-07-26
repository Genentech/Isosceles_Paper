import pysam

cell_barcode = "AAACCCAAGTATTGCC"

for sample_id in ["LIB5432309_SAM24385452", "LIB5432310_SAM24385453",
"LIB5432311_SAM24385454", "LIB5432312_SAM24385455", "LIB5432313_SAM24385456",
"LIB5432314_SAM24385457", "LIB5432315_SAM24385458", "LIB5432316_SAM24385459",
"LIB5427896_SAM24376275", "LIB5427897_SAM24376276"]:
    umi_iterator = itertools.product("ACGT", repeat=12)
    samfile_path = "bam/{0}.bam".format(sample_id)
    outfile_path = "bam/{0}_sicelore.bam".format(sample_id)
    with pysam.AlignmentFile(samfile_path, "rb") as samfile:
        with pysam.AlignmentFile(outfile_path, "wb", template=samfile) as outfile:
            for read in samfile:
                umi_sequence = "".join(next(umi_iterator))
                read.set_tag("BC", cell_barcode)
                read.set_tag("U8", umi_sequence)
                outfile.write(read)
