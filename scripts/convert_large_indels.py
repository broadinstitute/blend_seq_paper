import sys
import pysam


def replace_large_indels(cigar, threshold=50):
	for i, c in enumerate(cigar):
		if (c[0] == 1 or c[0] == 2) and c[1] > threshold:
			cigar[i] = (3, c[1])
	return cigar

with pysam.AlignmentFile(sys.argv[1], 'r') as bam:
	header = bam.header
	with pysam.AlignmentFile(sys.argv[2], 'wb', header=header) as out_bam:
		for read in bam:
			read.cigar = replace_large_indels(read.cigar)
			out_bam.write(read)
