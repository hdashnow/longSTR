# Parse TRF ngs dat file into tab delimited file
# Assumes TRF was run on aligned contigs.
# Each @ is a contig labeled with it's mapping position,
# each line is an STR

from strtools import normalise_str

trf_files = ['HG00733.h0.fasta.trf.dat',
    'HG00733.h1.fasta.trf.dat']

slop = 0 # Add this much slop to both sides of the variant position

def write_variant(outfile, this_variant):
    outfile.write(
        '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            this_variant['chrom'],
            this_variant['contig_start'] - slop, this_variant['contig_start'] + slop,
            this_variant['repeatunit'], this_variant['period'],
            this_variant['length_ru'], this_variant['length_bp'])
   )

for trf_file in trf_files:
    outfilename = trf_file + '.str'
    with open(outfilename, 'w') as outfile:
        #outfile.write('\t'.join(['chrom', 'contig_start', 'ref_stop', 'repeatunit', 'period', 'length_ru', 'length_bp']) + '\n')
        with open(trf_file) as trf_dat:
            for line in trf_dat:
                if line.startswith('@'):

                    # set new variant
                    splitlocus = line.strip('@').split(':')
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

                    this_variant = {'chrom': chrom, 'contig_start': contig_start, 'repeatunit': repeatunit,
                        'period': period, 'length_ru': length_ru, 'length_bp': length_bp}
                    write_variant(outfile, this_variant)

                    # consensus_size = splitline[4]
                    # percent_match = splitline[5]
                    # percent_indels = splitline[6]
                    # percent_A = splitline[8]
                    # percent_C = splitline[9]
                    # percent_G = splitline[10]
                    # percent_T = splitline[11]
                    # entropy = splitline[12]
                    # sequence = splitline[14]
                    # left_flank = splitline[15]
                    # right_flank = splitline[16]
