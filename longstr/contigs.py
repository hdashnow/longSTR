# Iterate over contigs in a bam file and calculate what base in the
# reference each base on the contig corresponds to

import pysam
from interlap import InterLap

# https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
# https://samtools.github.io/hts-specs/SAMv1.pdf
#               consumes_qry    consumes_ref
MATCH  = 0  # M yes             yes
INS    = 1  # I yes             no
DEL    = 2  # D no              yes
SKIP   = 3  # N no              yes
SOFT   = 4  # S yes             no
HARD   = 5  # H no              no
PAD    = 6  # P no              no
EQUAL  = 7  # = yes             yes
DIFF   = 8  # X yes             yes

def consumes_query(cigartuple):
    """If cigar operation consumes query return True, else return False"""
    if cigartuple[0] in [MATCH, INS, SOFT, EQUAL, DIFF]:
        return True
    else:
        return False

def consumes_reference(cigartuple):
    """If cigar operation consumes reference return True, else return False"""
    if cigartuple[0] in [MATCH, DEL, SKIP, EQUAL, DIFF]:
        return True
    else:
        return False

def get_intervals(sam):
    """Add all intervals from all contigs to an interval tree"""
    contig_intervals = InterLap()

    for contig in sam.fetch():
        qry_pos = 0
        ref_pos = contig.reference_start

        for cigartuple in contig.cigartuples:
            # current interval's query and ref start coords
            qry_start = qry_pos
            qry_end = qry_pos
            ref_start = ref_pos
            ref_end = ref_pos

            op_len = cigartuple[1]
            if consumes_query(cigartuple):
                qry_pos += op_len
                qry_end += op_len
            if consumes_reference(cigartuple):
                ref_pos += op_len
                ref_end += op_len

             # add the interval to the contig interval tree
            contig_intervals.add((qry_start, qry_end, (ref_start, ref_end)))
    return(contig_intervals)

def main():
    cramfile = "HG00512.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram"
    samfile = "/uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/git/STRling/working/chaisson_2019/data/HG00733.h0.chr1-1436179-1516179.sam"
    ref_fasta = ""
    trf_dat = "/uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/git/STRling/working/chaisson_2019/data/HG00733.h0.chr1-1436179-1516179.fasta.trf.dat"

    #sam = pysam.AlignmentFile(cramfile, "rc", )
    sam = pysam.AlignmentFile(samfile, "r")

    contig_intervals = get_intervals(sam)

    for i in contig_intervals:
        print(i)

if __name__ == "__main__":
    main()
