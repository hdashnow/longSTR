import os
import collections
import pysam
from interlap import InterLap

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
    return cigartuple[0] in [MATCH, INS, SOFT, EQUAL, DIFF]

def consumes_reference(cigartuple):
    """If cigar operation consumes reference return True, else return False"""
    return cigartuple[0] in [MATCH, DEL, SKIP, EQUAL, DIFF]


class Contig2Reference(object):
    def __init__(self, sam, reference=None):
        self.inter = collections.defaultdict(InterLap)
        self.add_sam(sam)
        self.reference = reference

    def add_sam(self, samfile):

        aligntype = os.path.splitext(samfile)[1]
        if aligntype == '.sam':
            sam = pysam.AlignmentFile(samfile, "r")
        elif aligntype == '.bam':
            sam = pysam.AlignmentFile(samfile, "rb")
        elif aligntype == '.cram':
            sam = pysam.AlignmentFile(samfile, "rc")
        else:
            sys.exit(f'File extension {aligntype} is not a recognized alignment format. Use .sam, .bam or .cram')

        for aln in sam:
            if aln.is_secondary or aln.is_unmapped or aln.mapping_quality < 40: continue 
            if aln.is_supplementary: continue # NOTE: not sure we shoulds skip these ?


            I = self.inter[aln.query_name]
            ref_name = aln.reference_name

            aln_off = 0
            ref_off = aln.reference_start
            for c in aln.cigartuples:
                op_len = c[1]
                cons_q = consumes_query(c)
                cons_r = consumes_reference(c)

                I.add((aln_off, aln_off + op_len * int(cons_q), (ref_name, ref_off, ref_off + op_len * int(cons_r))))

                aln_off += int(cons_q) * op_len
                ref_off += int(cons_r) * op_len

    def translate(self, read_name, pos):
        """translate a contig coordinate to a genomic coordinate. since there
        could be many mappings this is an iterator."""
        for (start, stop, (ref_chrom, ref_start, ref_stop)) in self.inter[read_name].find((pos, pos)):
            over = pos - start # how far past start
            #assert stop >= pos
            if stop < pos:
                print('read: ', read)
                print('stop: ', stop)
                print('pos: ', pos)
                sys.exit()
            yield (ref_chrom, ref_start + over)


if __name__ == "__main__":
    import sys
    c2r = Contig2Reference(sys.argv[1])

    for mapping in c2r.translate("000093F", 6862262):
        print("start:", mapping)
    for mapping in c2r.translate("000093F", 6862455):
        print("stop:", mapping)
