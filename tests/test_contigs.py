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

def test_translate():
    sam = "FRAXE.hap1.sam"

    c2r = Contig2Reference(sam)
    for mapping in c2r.translate("000093F", 6862262):
        print("start:", mapping)
        assert mapping == ('chrX', 148500605, False)
        break
    print()
    for mapping in c2r.translate("000093F", 6862455):
        print("stop:", mapping)
        assert mapping == ('chrX', 148500753, False)
        break

#def test_time_Contig2Reference():
#    import time
#    start = time.time()
#    c2r = Contig2Reference('/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/Porubsky_2019/HG00733.PUR.20191122_v1-1.HiFi.pg-rac2x.hap1.sam')
#    end = time.time()
#    print(end - start)
#    assert False
