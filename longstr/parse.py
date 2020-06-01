# Parse TRF ngs dat file into tab delimited file
# Assumes TRF was run on aligned contigs.
# Each @ is a contig labeled with it's mapping position,
# each line is an STR

import os
import sys
from .strtools import normalise_str

slop = 0 # Add this much slop to both sides of the variant position

def write_variant(outfile, variant):
    outfile.write(
        '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            variant['chrom'],
            variant['contig_start'] - slop, variant['contig_start'] + slop,
            variant['repeatunit'], variant['period'],
            variant['length_ru'], variant['length_bp'])
   )

def parse_dat(trf_file):
    with open(trf_file) as trf_dat:
        for line in trf_dat:
            if len(line.strip()) == 0: #XXX skip blank lines
                continue
            if line.startswith('@'):
                contig_id = line.strip('@')
                # set new variant
                splitlocus = contig_id.split(':')
                chrom = splitlocus[0]
                try:
                    splitlocus2 = splitlocus[1].split('/')
                    splitlocus3 = splitlocus2[0].split('-')
                    contig_start = int(splitlocus3[0])
                    contig_end = splitlocus3[1]
                except IndexError:
                    contig_start = 0
                    contig_end = 0
            else:
                splitline = line.split()
                start = splitline[0]
                end = splitline[1]
                period = splitline[2]
                length_ru = float(splitline[3])
                alignment_score = splitline[7]
                repeatunit = normalise_str(splitline[13])
                length_bp = int(end) - int(start) + 1

                #XXX adjust STR position for any upstream indels in the cotig? Track cigar string

                variant = {'chrom': chrom, 'contig_start': contig_start,
                    'repeatunit': repeatunit, 'period': period,
                    'length_ru': length_ru, 'length_bp': length_bp}
                yield (contig_id, variant)

def parse_trf(trf_file, outfilename = ''):
    if outfilename != '':
        outfile = open(outfilename, 'w')

    for contig_id, variant in parse_dat(trf_file):
        # Adjust position to reference coordinates
        if outfilename != '':
            write_variant(outfile, variant)
