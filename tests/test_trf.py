import sys
sys.path.append("../longstr/")
from trf import *
import pytest
from interlap import InterLap

def test_parse_dat():
    trf_contigs = parse_dat('test.trf.dat')
    i = 0
    for contig, variants in trf_contigs:
        i += 1
        if i == 1:
            assert contig == 'chr1:1000-1350/0'
            #assert list(variant_intervals.find((95, 100)))[0][0:2] == (91, 139)
        if i == 2:
            assert contig == 'chr2:1000-1350/0'
        print(contig)
        for variant in variants:
            print(variant)

@pytest.mark.parametrize("pos, expected", [
    (25, 1024),
    (55, 1049),
    (50, 1049),
    (150, None),
])
def test_get_ref_pos(pos, expected):
    intervals = InterLap()    
    intervals.add( (0, 50, (999, 1049)) )
    intervals.add( (50, 100, (1049, 1049)) )
    assert get_ref_pos(pos, intervals) == expected

def test_trf_to_genome():
    trf_to_genome('test.sam', 'test.trf.dat', 'test_output.txt')
