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

def test_consumes_reference():
    sam = pysam.AlignmentFile("test.sam", "r")
    contig_intervals = get_intervals(sam)
    for i in contig_intervals:
        assert i == (0, 50, (999, 1049))
        break
