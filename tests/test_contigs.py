import sys
sys.path.append("../longstr/")
from contigs import *
import pytest

@pytest.mark.parametrize("cigartuple, expected", [
    ((0, 10), True),
    ((2, 10), False),
])
def test_consumes_query(cigartuple, expected):
    assert consumes_query(cigartuple) == expected

def test_get_intervals():
    sam = pysam.AlignmentFile("test.sam", "r")
    # Cigar: 50=50I100=100D100=
    # Aligned to chr1 1000
    for contig_id, location, contig_intervals in get_intervals(sam):
        print(contig_id)
        for i in contig_intervals:
            print(i)
        if contig_id == 'chr2:1000-1350/0':
            assert list(contig_intervals.find((0, 21))) == [(0, 50, (999, 1049))]
            assert list(contig_intervals.find((55, 70))) == [(50, 100, (1049, 1049))]
            assert list(contig_intervals.find((90, 110))) == [(50, 100, (1049, 1049)),
                                                                (100, 200, (1049, 1149))]
        
        # May need to think about this one:
        #assert list(contig_intervals.find((200, 200))) == [(200, 300, (1249, 1349))]

# Check supplementary alignments are correctly included
def test_get_intervals_supp():
    sam = pysam.AlignmentFile("test.sam", "r")
    # Cigar: 50=50I100=100D100=
    # Aligned to chr1 1000 and chr2 3000
    for contig_id, location, contig_intervals in get_intervals(sam):
        print(contig_id)
        for i in contig_intervals:
            print(i)
            print()
            for j in contig_intervals.find((0, 21)):
                print(j)
        assert list(contig_intervals.find((0, 21))) == [(0, 50, (2999, 3049)),
                                                        (0, 50, (999, 1049))]
        break
