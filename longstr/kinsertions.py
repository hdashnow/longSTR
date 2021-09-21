#!/usr/bin/env python
"""Find insertions containing STR sequence using k-mer counting
"""

import argparse
from cyvcf2 import VCF
import nimporter, kmers
import sys

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('VCF', type=str)
    parser.add_argument('OUT', type=str, default=sys.stdout,
                            help='Output tsv filename')
    parser.add_argument('--min_allele', type=int, default=10,
                            help='Minimum insertion size to check (default: %(default)s)',)


    return parser.parse_args()

def count_kmers(s):
    freq_kmer = 'GT'
    kmer_count = s.count('GT')

    # make a new ktable
    ksize = 2
    target_table_size = 5e8
    num_tables = 4
    ktable = khmer.Counttable(ksize, target_table_size, num_tables)

    # count all k-mers in the given string
    ktable.consume(s)

    # run through all entries. if they have nonzero presence, print.
    for kmer in ktable.get_kmers(s):
        print(kmer)
    return (freq_kmer, kmer_count)
 
def report_kmers(vcffile, outname, min_allele_len = 10):

    with open(outname, 'w') as outfile:
        header = ['chrom', 'start', 'end', 
                    'repeatunit', 'prop_repeat', 'ref', 'alt']
        outfile.write('\t'.join(header) + '\n')
        for variant in VCF(vcffile):
            for allele in variant.ALT:
                if len(allele) < min_allele_len:
                    continue
                # Get most frequest k-mer and the proption of the allele it covers
                kmer, prop = kmers.kmer_freq(allele)
                # Report genomic coordinates of STR insertions
                out_list = [variant.CHROM, variant.start, variant.end, 
                            kmer, prop, variant.REF, allele]
                out_list = [str(x) for x in out_list]
                outfile.write('\t'.join(out_list) + '\n')

def main():
    args = parse_args()

    vcffile = '/uufs/chpc.utah.edu/common/HIPAA/u6026198/storage/git/longSTR/working/test.vcf'
    outname = 'test-out.tsv'

    report_kmers(args.VCF, args.OUT, args.min_allele)

if __name__ == '__main__':
    main()
