# this script is adapted from https://github.com/vierstralab/WASP/blob/master/mapping/rmdup.py

import pysam
import argparse
from numpy.random import default_rng


rng = default_rng(42)

parser = argparse.ArgumentParser()
parser.add_argument('input_bam', help="input BAM or SAM file (must be sorted!)")
parser.add_argument("output_bam", help="output BAM or SAM file")

options = parser.parse_args()

if options.input_bam.endswith(".sam") or options.input_bam.endswith("sam.gz"):
    infile = pysam.Samfile(options.input_bam, "r")
else:
    # assume binary BAM file
    infile = pysam.Samfile(options.input_bam, "rb")

if options.output_bam.endswith(".sam"):
    # output in text SAM format
    outfile = pysam.Samfile(options.output_bam, "w", template=infile)
elif options.output_bam.endswith(".bam"):
    # output in binary compressed BAM format
    outfile = pysam.Samfile(options.output_bam, "wb", template=infile)
else:
    raise ValueError("name of output file must end with .bam or .sam")


chr = pos = None
read_cache = []
flags = []
for read in infile:
    print(chr, pos)
    if read.rname != chr or read.pos != pos:
        if len(read_cache) > 0:
            outfile.write(rng.choice(read_cache))
        if len(flags) > 1:
            print(flags)
        chr = read.rname
        pos = read.pos
        read_cache = []
        flags = []
    read_cache.append(read)
    flags.append(read.flag)
infile.close()
outfile.close()
