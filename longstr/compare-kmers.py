# Compare STRling calls to PacBio insertion calls with k-mer content

import argparse
import sys
import os
from strtools import normalise_str
import pandas as pd
import pyranges as pr

def parse_args():
    # top-level parser
    parser = argparse.ArgumentParser()

    parser.add_argument('--strling', type=str,
                        help='STRling outlier tsv file for one sample')
    parser.add_argument('--pacbio', type=str, nargs='+',
                        help='STR variants from PacBio with kmers in bed format (two, one for each haplotype or primary and associate contigs)')
    parser.add_argument('--cov', type=str, nargs='+',
                        help='Contig coverage in bed format (two, one for each haplotype or primary and associate contigs)')
    parser.add_argument('--centromeres', type=str, help='Optional bed file of centromeres to annotate')
    parser.add_argument('--telomeres', type=str, help='Optional bed file of telomeres to anotate')
    parser.add_argument('--LCRs', type=str, help='Optional bed file of Low copy repeats to anotate')
    parser.add_argument('--slop', type=int, default=500,
                        help='Consider variants the same if they are this distance apart (and have the same repeat unit)')
    parser.add_argument('--out', type=str, default = '',
                        help='output file of variants in annotated bed format')

    return parser.parse_args()

def get_sample(fullpath):
    """Get the sample ID from the filename"""
    basename = os.path.basename(fullpath)
    return(basename.rsplit('-', maxsplit = 1)[0].split('.')[0])

def parse_bed(filename):
    """Parse bed file with header row starting with #"""
    sample_id = get_sample(filename)
    try:
        df = pd.read_csv(filename, sep = '\t')
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    df.columns = [col.strip('#') for col in df.columns]

    for index, row in df.iterrows():
        repeatunit = row['repeatunit']
        if pd.isna(repeatunit):
            sys.stderr.write(str(row)+'\n')
            sys.exit('ERROR: NA value found in the repeatunit column of {}: {}\n'.format(filename, repeatunit))
        for base in repeatunit:
            if base not in ['A', 'T', 'C', 'G']:
                sys.stderr.write(str(row)+'\n')
                sys.exit('ERROR: Non-DNA found in the repeatunit column of {}: {}\n'.format(filename, repeatunit))
    df['sample'] = sample_id

    # Rename columns to match pyranges expectations
    df = df.rename(columns={'chrom': 'Chromosome',
                    'left': 'Start', 'ref_start': 'Start',
                    'start': 'Start', 'end': 'End',
                    'right': 'End', 'ref_end': 'End'})

    return(df)

def prop_str(S):
    """Takes a pd Series with at least the indices 'alt' and 'repeatunit',
        both strings. Return the number of occurances of repeatunit in alt"""

    if S['alt'] is None:
        return 0

    count = S['alt'].count(S['repeatunit'])

    return count * len(S['repeatunit']) / len(S['alt'])

def count_str(S):
    """Takes a pd Series with at least the indices 'alt' and 'repeatunit',
        both strings. Return the number of occurances of repeatunit in alt"""

    if S['alt'] is None:
        return 0

    count = S['alt'].count(S['repeatunit'])

    return count

def count_indel(S):
    """Takes a pd Series with at least the indices 'alt' and 'ref',
        both strings. Return len(alt) - len(ref)"""

    if S['alt'] is None:
        return 0

    indel = len(S['alt']) - len(S['ref'])

    return indel


def match_closest(strling_df, pacbio_hap1_df, pacbio_hap2_df, slop=500):
    """Annotate with the closest locus"""

    if strling_df.empty:
        return strling_df

    for i, pacbio_df in zip((1,2), (pacbio_hap1_df, pacbio_hap2_df)):

        if pacbio_df.empty:
            continue #XXX May need to add extra columns?
    
        strling_pr = pr.PyRanges(strling_df)
        pacbio_pr = pr.PyRanges(pacbio_df)
    
        # Annotate with the closest locus
        nearest_pr = strling_pr.nearest(pacbio_pr)
        nearest_df = nearest_pr.df

        if nearest_df.empty:
            continue

        # Remove pacbio variants more than slop bp away
        nearest_columns = ['Start_b', 'End_b', 'repeatunit_norm_b', 'period', 'length_ru', 
                            'length_bp', 'prop_repeat', 'sample_b', 'Distance', 'ref', 'alt']
        nearest_df.loc[nearest_df.Distance > slop, nearest_columns] = None

        # Create a unique index
        strling_df['locus'] = strling_df['Chromosome'] + '-' + strling_df['Start'].astype(str
            ) + '-' + strling_df['End'].astype(str) + '-' + strling_df['repeatunit']
        strling_df.set_index('locus', inplace = True)
        nearest_df['locus'] = nearest_df['Chromosome'].astype(str) + '-' + nearest_df['Start'
            ].astype(str) + '-' + nearest_df['End'].astype(str) + '-' + nearest_df['repeatunit']
        nearest_df.set_index('locus', inplace = True)
        nearest_df['strlingRUcount'] = nearest_df.apply(count_str, axis=1) #rows
        nearest_df['strlingRUprop'] = nearest_df.apply(prop_str, axis=1) #rows
        nearest_df['indel'] = nearest_df.apply(count_indel, axis=1) #rows

        nearest_df = nearest_df.filter(['repeatunit_norm_b', 'prop_repeat', 
                                'Distance', 'strlingRUprop', 'strlingRUcount', 'indel'])
        nearest_df.rename(columns={'repeatunit_norm_b': f'repeatunit_hap{i}',
                                    'prop_repeat': f'prop_repeat_hap{i}',
                                    'strlingRUprop': f'strlingRUprop_hap{i}',
                                    'strlingRUcount': f'strlingRUcount_hap{i}',
                                    'Distance': f'Distance_hap{i}',
                                    'indel': f'indel_hap{i}',
                                    }, inplace = True)
        strling_df = strling_df.merge(nearest_df, how = 'left', left_index = True,
                                        right_index = True)

    return strling_df

def match_variants(strling_df, pacbio_hap1_df, pacbio_hap2_df, slop):
    """Match up pacbio pacbio variants to their corresponding strling variants
        - Within X bp of slop
        - Same repeat unit
        Any leftover pacbio pacbio variants added on to the end"""

    strling_df['repeatunit_norm'] = strling_df['repeatunit'].apply(normalise_str)
    pacbio_hap1_df['repeatunit_norm'] = pacbio_hap1_df['repeatunit'].apply(normalise_str)
    pacbio_hap2_df['repeatunit_norm'] = pacbio_hap2_df['repeatunit'].apply(normalise_str)

    # Annotate with closest insertion (any repeatunit)
    all_closest_df = match_closest(strling_df, pacbio_hap1_df, pacbio_hap2_df, slop)
    
    return all_closest_df

def annotate_cov(strling_df, cov_hap1_pr, cov_hap2_pr):
    """Annotate strling calls with contig coverage overlap"""
    strling_pr = pr.PyRanges(strling_df)

    for i, hap_pr in zip((1,2), (cov_hap1_pr, cov_hap2_pr)):
        cov_pr = strling_pr.coverage(hap_pr)
        cov_df = cov_pr.df

        # If FractionOverlaps < 1, set NumberOverlaps to 0 
        cov_df.loc[cov_df.FractionOverlaps < 1, 'NumberOverlaps'] = 0

       # Create a unique index
        cov_df['locus'] = cov_df['Chromosome'].astype(str) + '-' + cov_df['Start'
            ].astype(str) + '-' + cov_df['End'].astype(str) + '-' + cov_df['repeatunit']
        cov_df.set_index('locus', inplace = True)

        cov_df = cov_df.filter(['NumberOverlaps'])
        cov_df.rename(columns={'NumberOverlaps': f'Hap{i}Cov'}, inplace = True)
        strling_df = strling_df.merge(cov_df, how = 'left', left_index = True,
                                        right_index = True)
        strling_df[f'Hap{i}Cov'] = strling_df[f'Hap{i}Cov'].fillna(0)

    return(strling_df)

def annotate_feature_cov(strling_df, bed, feature):
    """Annotate strling calls with if they overlap a feature"""
    strling_pr = pr.PyRanges(strling_df)
    bed_pr = pr.readers.read_bed(bed)

    cov_pr = strling_pr.coverage(bed_pr)
    cov_df = cov_pr.df

    # Create a unique index
    cov_df['locus'] = cov_df['Chromosome'].astype(str) + '-' + cov_df['Start'
        ].astype(str) + '-' + cov_df['End'].astype(str) + '-' + cov_df['repeatunit']
    cov_df.set_index('locus', inplace = True)

    cov_df = cov_df.filter(['NumberOverlaps'])
    cov_df.rename(columns={'NumberOverlaps': f'{feature}Cov'}, inplace = True)
    strling_df = strling_df.merge(cov_df, how = 'left', left_index = True,
                                    right_index = True)
    strling_df[f'{feature}Cov'] = strling_df[f'{feature}Cov'].fillna(0)

    return(strling_df)

def main():
    args = parse_args()
    if len(args.pacbio) != 2:
        exit(f'ERROR: Expected 2 pacbio files, got {len(args.pacbio)}: {args.pacbio}')
    if len(args.cov) != 2:
        exit(f'ERROR: Expected 2 cov files, got {len(args.cov)}: {args.cov}')

    # Parse inputs as pandas data frame (df) or pyranges (pr) objects
    strling_df = pd.read_csv(args.strling, sep='\t')
    strling_df.columns = [col.strip('#') for col in strling_df.columns]
    strling_df = strling_df.rename(columns={'chrom': 'Chromosome',
                    'left': 'Start', 'right': 'End'})
    # +1 end of any regions of length 0 (Start == End) to avoid pyranges error
    strling_df.loc[strling_df['Start'] == strling_df['End'], 'End'] += 1


    pacbio_hap1_df = parse_bed(args.pacbio[0])
    pacbio_hap2_df = parse_bed(args.pacbio[1])

    cov_hap1_pr = pr.read_bed(args.cov[0])
    cov_hap2_pr = pr.read_bed(args.cov[1])

    # Annotate strling calls with corresponding PacBio calls
    strling_df = match_variants(strling_df, pacbio_hap1_df, pacbio_hap2_df, args.slop)

    # Annotate strling calls with coverage overlap
    strling_df = annotate_cov(strling_df, cov_hap1_pr, cov_hap2_pr)

    # Annotate with additional stuff if requested
    if args.centromeres:
        strling_df = annotate_feature_cov(strling_df, args.centromeres, 'centromere')
    if args.telomeres:
        strling_df = annotate_feature_cov(strling_df, args.telomeres, 'telomere')
    if args.LCRs:
        strling_df = annotate_feature_cov(strling_df, args.LCRs, 'LCR')

    strling_df.to_csv(args.out, sep='\t', index=False)

if __name__ == '__main__':
    main()
