# Parse TRF ngs dat file into tab delimited file
# Assumes TRF was run on aligned contigs.
# Each @ is a contig labeled with it's mapping position,
# each line is an STR

import argparse
import os
import sys
from strtools import normalise_str
from contigs import Contig2Reference
import itertools
import pysam
import statistics as stats

def parse_args():
    # top-level parser
    parser = argparse.ArgumentParser()

    # subcommands
    parser.add_argument('BAM', type=str,
                        help='bam or cram file')
    #parser.add_argument('--fasta', type=str,
    #                    help='reference genome fasta (required for cram)')
    parser.add_argument('--dat', type=str,
                        help='TRF dat file')
    parser.add_argument('--out', type=str, default = '',
                        help='output file of variants in annotated bed format')

    return parser.parse_args()

def parse_dat(trf_file):
    with open(trf_file) as trf_dat:
        for line in trf_dat:
            if len(line.strip()) == 0: # skip blank lines
                continue
            if line.startswith('@'):
                # Complete previous variant
                try:
                    yield (contig_id, variants)
                except UnboundLocalError: # on first line, since not yet populated
                    pass

                contig_id = line.strip('@').strip().split()[0]
                variants = []
                # set new variant
            else:
                splitline = line.split()
                start = int(splitline[0])
                end = int(splitline[1])
                period = splitline[2]
                length_ru = float(splitline[3])
                alignment_score = splitline[7]
                repeatunit = normalise_str(splitline[13])
                length_bp = end - start + 1

                variant = {'start': start, 'end': end,
                            'repeatunit': repeatunit, 'period': period,
                            'length_ru': length_ru, 'length_bp': length_bp}
                variants.append(variant)
        # Complete final variant
        yield (contig_id, variants)

def get_ref_pos(pos, intervals):
    results = list(intervals.find((pos, pos)))
    if len(results) == 0:
        return
    all_results = []
    for result in results:
        result_contig = (result[0], result[1])
        result_ref = result[2]
        # Result reference is 1 bp
        if result_ref[0] == result_ref[1]:
            all_results.append(result_ref[0])
        # Result is a range so need to calculate position in range
        else:
            # If both ranges have a non-zero lenth, check the are equal
            if result_contig[0] != result_contig[1]:
                try:
                    assert result_ref[1] - result_ref[0] == result_contig[1] - result_contig[0]
                except AssertionError:
                    sys.stderr.write('WARNING, inconsistent range lengths:\n')
                    sys.stderr.write(f'result_ref {result_ref} len: {result_ref[1] - result_ref[0]}\n')
                    sys.stderr.write(f'result_contig {result_contig} len: {result_contig[1] - result_contig[0]}\n')
                    sys.stderr.write(f'result would be {result_ref[0] + pos - result_contig[0]}\n')
                    assert False
            all_results.append(result_ref[0] + pos - result_contig[0])
    # Check all results are equal
    if len(set(all_results)) == 1:
        return all_results[0]
    else:
        try:
            return stats.mode(all_results)
        except stats.StatisticsError:
            return round(stats.mean(all_results))

header = '\t'.join([
            '#chrom',
            'ref_start', 'ref_end',
            'repeatunit', 'period',
            'length_ru', 'length_bp',
            'indel',
            ]) + '\n'

def write_variant(outfile, variant):
    outfile.write(
        '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            variant['chrom'],
            variant['ref_start'], variant['ref_end'],
            variant['repeatunit'], variant['period'],
            variant['length_ru'], variant['length_bp'],
            variant['indel'],
        )
   )

def get_ref_pos(contigID, position, contig2ref):
    chrom_matches = []
    pos_matches = []
    for i, mapping in enumerate(contig2ref.translate(contigID, position)):
        chrom_matches.append(mapping[0])
        pos_matches.append(mapping[1])
    if len(chrom_matches) == 0:
        return (None, None)
    # Check all results are equal and if not chose one
    assert len(set(chrom_matches)) == 1
    #    print(chrom_matches)
    #    print(pos_matches)
    if len(set(pos_matches)) == 1:
        return (chrom_matches[0], pos_matches[0])
    else:
        try:
            return (chrom_matches[0], stats.mode(pos_matches))
        except stats.StatisticsError:
            return (chrom_matches[0], round(stats.mean(pos_matches)))

def trf_to_genome(samfile, trf_file, outfilename = ''):
    """Read in trf calls on contigs and the corresponding alignments and use these
    to transfer the calls from contig coordinates to reference genome coordinates.
    """
    if outfilename != '':
        outfile = open(outfilename, 'w')
        outfile.write(header)

    contig2ref = Contig2Reference(samfile)

    for contigID, variants in parse_dat(trf_file):

        for variant in variants:
            chrom_start, variant['ref_start'] = get_ref_pos(contigID, variant['start'], contig2ref)
            chrom_end, variant['ref_end'] = get_ref_pos(contigID, variant['end'], contig2ref)
            if chrom_start == None or chrom_start == None: # Skip variants with no match in the sam
                continue
            assert chrom_start == chrom_end
            #if chrom_start != chrom_end:
            #    print(variant)
            variant['chrom'] = chrom_start
            variant['indel'] = (variant['start'] - variant['end']) - (variant['ref_start'] - variant['ref_end'])
            if outfilename != '':
                write_variant(outfile, variant)

def main():
    args = parse_args()

    trf_to_genome(args.BAM, args.dat, args.out)

if __name__ == '__main__':
    main()
