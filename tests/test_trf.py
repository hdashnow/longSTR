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

def test_trf_to_genome():
    trf_to_genome('test.sam', 'test.trf.dat', 'test_output.txt')
