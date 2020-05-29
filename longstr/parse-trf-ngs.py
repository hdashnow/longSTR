from strtools import normalise_str

trf_files = ['HG00514.insertions.fasta.ngs.dat',
    'HG00733.insertions.fasta.ngs.dat', 'NA19240.insertions.fasta.ngs.dat']

slop = 100 # Add this much slop to both sides of the variant position

for trf_file in trf_files:
    outfilename = trf_file + '.str'
    with open(outfilename, 'w') as outfile:
        #outfile.write('\t'.join(['chrom', 'ref_start', 'ref_stop', 'repeatunit', 'period', 'length_ru', 'length_bp', 'total_ins_size']) + '\n')
        with open(trf_file) as trf_dat:
            this_variant = {} #necessary?
            for line in trf_dat:
                if line.startswith('@'):

                    # write out previous variant
                    for repeatunit in this_variant:
                        outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(this_variant[repeatunit]['chrom'],
                        this_variant[repeatunit]['ref_start'] - slop, this_variant[repeatunit]['ref_start'] + slop,
                        repeatunit, this_variant[repeatunit]['period'],
                        this_variant[repeatunit]['length_ru'], this_variant[repeatunit]['length_bp'],
                        this_variant[repeatunit]['total_ins_size']))
                    # set new variant
                    this_variant = {}
                    splitlocus = line.strip('@').split(':')
                    chrom = splitlocus[0]
                    splitlocus2 = splitlocus[1].split('|')
                    ref_start = int(splitlocus2[0])
                    total_ins_size = splitlocus2[1]
                else:
                    splitline = line.split()
                    start = splitline[0]
                    end = splitline[1]
                    period = splitline[2]
                    length_ru = float(splitline[3])
                    alignment_score = splitline[7]
                    repeatunit = normalise_str(splitline[13])
                    length_bp = int(end) - int(start) + 1

                    if repeatunit not in this_variant:
                        this_variant[repeatunit] = {'chrom': chrom, 'ref_start': ref_start,
                        'period': period, 'length_ru': length_ru, 'length_bp': length_bp, 'total_ins_size': total_ins_size}
                    else:
                        this_variant[repeatunit]['length_ru'] += length_ru
                        this_variant[repeatunit]['length_bp'] += length_bp
            # write out final variant
            for repeatunit in this_variant:
                outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(this_variant[repeatunit]['chrom'],
                this_variant[repeatunit]['ref_start'] - slop, this_variant[repeatunit]['ref_start'] + slop,
                repeatunit, this_variant[repeatunit]['period'],
                this_variant[repeatunit]['length_ru'], this_variant[repeatunit]['length_bp'],
                this_variant[repeatunit]['total_ins_size']))

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
