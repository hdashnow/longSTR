import argparse
import sys
import os
from longstr.parse import parse_trf

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def run_trf(args):
    print(args)

def run_parse(args):
    parse_trf(args.dat, args.out)

def parse_args():
    # top-level parser
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title='subcommands')   
 
    # subcommands
    parser_trf = subparsers.add_parser('trf', help='Run Tandem Repeats Finder on bam/cram')
    parser_trf.add_argument('BAM', type=str,
                            help='bam or cram file')
    parser_trf.add_argument('--fasta', type=str,
                            help='reference genome fasta (required for cram)')
    parser_trf.set_defaults(func=run_trf)
    
    parser_parse = subparsers.add_parser('parse', help='parse TRF output')
    parser_parse.add_argument('BAM', type=str,
                            help='bam or cram file')
    parser_parse.add_argument('--fasta', type=str,
                            help='reference genome fasta (required for cram)')
    parser_parse.add_argument('--dat', type=str,
                            help='TRF dat file (inferred from BAM name by default)')
    parser_parse.add_argument('--out', type=str, default = '',
                            help='output file of variants in annotated bed format')
    parser_parse.set_defaults(func=run_parse)

    return parser.parse_args()

def main():
    args = parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
