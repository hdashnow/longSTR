# Parse TRF ngs dat file into tab delimited file
# Assumes TRF was run on aligned contigs.
# Each @ is a contig labeled with it's mapping position,
# each line is an STR

import os
import sys
from strtools import normalise_str
from contigs import get_intervals
import itertools
import pysam

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

                contig_id = line.strip('@').strip()
                variants = []
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
                start = int(splitline[0])
                end = int(splitline[1])
                period = splitline[2]
                length_ru = float(splitline[3])
                alignment_score = splitline[7]
                repeatunit = normalise_str(splitline[13])
                length_bp = end - start + 1

                variant = {'chrom': chrom, 'start': start, 'end': end,
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
            assert result_ref[1] - result_ref[0] == result_contig[1] - result_contig[0]
            all_results.append(result_ref[0] + pos - result_contig[0])
    # Check all results are equal
    if len(set(all_results)) == 1:
        return all_results[0]
    else:
        sys.exit('Inconsistent results for position, error?')

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

def trf_to_genome(alignfile, trf_file, outfilename = ''):
    """Read in trf calls on contigs and the corresponding alignments and use these
    to transfer the calls from contig coordinates to reference genome coordinates.
    """
    if outfilename != '':
        outfile = open(outfilename, 'w')
        outfile.write(header)

    aligntype = os.path.splitext(alignfile)[1]
    if aligntype == '.sam':
        sam = pysam.AlignmentFile(alignfile, "r")
    elif aligntype == '.bam':
        sam = pysam.AlignmentFile(alignfile, "rb")
    elif aligntype == '.cram':
        sam = pysam.AlignmentFile(alignfile, "rc")
    else:
        sys.exit(f'File extension {aligntype} is not a recognized alignment format. Use .sam, .bam or .cram')

    contig_data = get_intervals(sam)
    variant_data = parse_dat(trf_file)
    # Iterate through both at once
    for contig_tuple, variant_tuple in itertools.zip_longest(contig_data,
                                        variant_data, fillvalue=None):
        if contig_tuple == None or variant_tuple == None:
            sys.exit('Contigs differ between bam and trf dat file. Check same was used for both')
        if contig_tuple[0] != variant_tuple[0]:
            sys.exit('ID mismatch: {contig_tuple[0]}, {variant_tuple[0]}')
        for interval in contig_tuple[1]:
            contig_intervals = contig_tuple[1]
        for variant in variant_tuple[1]:
            variant['ref_start'] = get_ref_pos(variant['start'], contig_intervals)
            variant['ref_end'] = get_ref_pos(variant['end'], contig_intervals)
            variant['indel'] = (variant['start'] - variant['end']) - (variant['ref_start'] - variant['ref_end'])
            if outfilename != '':
                write_variant(outfile, variant)

