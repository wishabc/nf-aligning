import argparse
from numpy.random import default_rng
import pysam
import shutil

rng = default_rng(42)

def main(infile, outfile, reads_count, sample_size):
    # if total number of reads less than number of reads to select,
    # then copy input file to output file
    if sample_size >= reads_count:
        shutil.copyfile(infile, outfile)
        return

    sorted_read_indexes = set(rng.choice(range(reads_count), sample_size, replace=False))

    print('Selecting %d read pairs' % len(sorted_read_indexes))

    with pysam.AlignmentFile(infile, 'rb') as in_alignment_file, pysam.AlignmentFile(outfile, 'wb', template=in_alignment_file) as out_alignment_file:
        for i, alignment in enumerate(in_alignment_file):
            if i in sorted_read_indexes:
                out_alignment_file.write(alignment)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("outfile")
    parser.add_argument("paired_reads_count")
    parser.add_argument("paired_reads_count_to_select")
    args = parser.parse_args()

    # bam file to read from
    infile = args.infile
    # bam file to write to
    outfile = args.outfile

    main(args.infile, args.outfile,
        int(args.paired_reads_count), int(args.paired_reads_count_to_select))