# Iterate over contigs in a bam file and calculate what base in the
# reference each base on the contig corresponds to

import pysam

"""
cigartuples             consumes reference
M   BAM_CMATCH  0       yes
I   BAM_CINS    1       no
D   BAM_CDEL    2       yes
N   BAM_CREF_SKIP   3   yes
S   BAM_CSOFT_CLIP  4   no
H   BAM_CHARD_CLIP  5   no
P   BAM_CPAD    6       no
=   BAM_CEQUAL  7       yes
X   BAM_CDIFF   8       yes
B   BAM_CBACK   9       ?
"""

def consumes_ref(cigartuple):
    """If cigar operation consumes reference return True, else return False"""
    if cigartuple[0] in [0, 2, 3, 7, 8]:
        return True
    else:
        return False

def main():
    cramfile = "HG00512.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram"
    ref_fasta = ""
    trf_dat = "/uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/git/STRling/working/chaisson_2019/data/HG00733.h0.chr1-1436179-1516179.fasta.trf.dat"

    #sam = pysam.AlignmentFile(cramfile, "rc", )
    sam = pysam.AlignmentFile("HG00733.h0.chr1-1436179-1516179.sam", "r")

    for read in sam.fetch():
        footprint = 0
        for cigartuple in read.cigartuples:
            if consumes_ref(cigartuple):
                footprint += cigartuple[1]
        print(footprint)
        print(read.query_alignment_length)
        exit()

if __name__ == "__main__":
    main()
