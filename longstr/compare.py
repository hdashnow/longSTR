# Parse TRF ngs dat file into tab delimited file
# Assumes TRF was run on aligned contigs.
# Each @ is a contig labeled with it's mapping position,
# each line is an STR

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
                        help='STRling outlier tsv file')
    parser.add_argument('--trf', type=str,
                        help='STR variants from TRF in bed format')
    parser.add_argument('--cov', type=str,
                        help='Contig coverage in bed format')
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
        df = pd.read_csv(filename, delim_whitespace = True)
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    df.columns = [col.strip('#') for col in df.columns]

    for repeatunit in df['repeatunit']:
        for base in repeatunit:
            if base not in ['A', 'T', 'C', 'G']:
                sys.exit('ERROR: Non-DNA found in the third column of {}: {}\n'.format(filename, repeatunit))
    df['sample'] = sample_id

    # Rename columns to match pyranges expectations
    df = df.rename(columns={'chrom': 'Chromosome',
                    'left': 'Start', 'ref_start': 'Start',
                    'right': 'End', 'ref_end': 'End'})

    return(df)

def match_closest(strling_df, trf_df, this_repeatunit, slop=200):
    """Filter to a specific repeat unit then annotate with the closest locus"""
    strling_df = strling_df.loc[strling_df['repeatunit_norm'] == this_repeatunit].copy()

    if strling_df.empty:
        return strling_df

    trf_df = trf_df[trf_df['repeatunit_norm'] == this_repeatunit]
    if trf_df.empty:
        return #XXX May need to add extra columns?

    strling_pr = pr.PyRanges(strling_df)
    trf_pr = pr.PyRanges(trf_df)

    # Annotate with the closest locus
    nearest_pr = strling_pr.nearest(trf_pr)
    nearest_df = nearest_pr.df

    if nearest_df.empty:
        return

    # Remove pacbio variants more than slop bp away
    nearest_columns = ['Start_b', 'End_b', 'repeatunit_norm_b', 'period', 'length_ru', 'length_bp', 'indel', 'sample_b', 'Distance']
    nearest_df.loc[nearest_df.Distance > slop, nearest_columns] = None

    # Create a unique index
    strling_df['locus'] = strling_df['Chromosome'] + '-' + strling_df['Start'].astype(str
        ) + '-' + strling_df['End'].astype(str) + '-' + strling_df['repeatunit']
    strling_df.set_index('locus', inplace = True)
    nearest_df['locus'] = nearest_df['Chromosome'].astype(str) + '-' + nearest_df['Start'
        ].astype(str) + '-' + nearest_df['End'].astype(str) + '-' + nearest_df['repeatunit']
    nearest_df.set_index('locus', inplace = True)

    nearest_df = nearest_df.filter(['repeatunit_norm_b', 'indel', 'Distance'])
    nearest_df.rename(columns={'repeatunit_norm_b': f'repeatunit_hap',
                                'indel': f'indel_hap',
                                'Distance': f'Distance_hap'}, inplace = True)
    strling_df = strling_df.merge(nearest_df, how = 'left', left_index = True,
                                    right_index = True)

    return strling_df

def match_variants(strling_df, trf_hap1_df, slop):
    """Match up pacbio trf variants to their corresponding strling variants
        - Within X bp of slop
        - Same repeat unit
        Any leftover pacbio trf variants added on to the end"""

    strling_df['repeatunit_norm'] = strling_df['repeatunit'].apply(normalise_str)
    trf_hap1_df['repeatunit_norm'] = trf_hap1_df['repeatunit'].apply(normalise_str)

    # Break down by repeat unit before comparing
    all_closest_df = pd.DataFrame()
    for this_ru in set(strling_df['repeatunit_norm']):
        all_closest_df = all_closest_df.append(match_closest(strling_df, trf_hap1_df, this_ru, slop))
    
    return all_closest_df

def annotate_cov(strling_df, hap_pr):
    """Annotate strling calls with contig coverage overlap"""
    strling_pr = pr.PyRanges(strling_df)

    cov_pr = strling_pr.coverage(hap_pr)
    cov_df = cov_pr.df

    # If FractionOverlaps < 1, set NumberOverlaps to 0
    cov_df.loc[cov_df.FractionOverlaps < 1, 'NumberOverlaps'] = 0

    # Create a unique index
    strling_df['locus'] = strling_df['Chromosome'] + '-' + strling_df['Start'].astype(str
        ) + '-' + strling_df['End'].astype(str) + '-' + strling_df['repeatunit']
    strling_df.set_index('locus', inplace = True)
    cov_df['locus'] = cov_df['Chromosome'].astype(str) + '-' + cov_df['Start'
        ].astype(str) + '-' + cov_df['End'].astype(str) + '-' + cov_df['repeatunit']
    cov_df.set_index('locus', inplace = True)

    cov_df = cov_df.filter(['NumberOverlaps'])

    cov_df.rename(columns={'NumberOverlaps': f'HapCov'}, inplace = True)
    strling_df = strling_df.merge(cov_df, how = 'left', left_index = True,
                                        right_index = True)
    return(strling_df)

def main():
    args = parse_args()

    # Parse inputs as pandas data frame (df) or pyranges (pr) objects
    strling_df = pd.read_csv(args.strling, sep='\t')
    strling_df = strling_df.rename(columns={'chrom': 'Chromosome',
                    'left': 'Start', 'right': 'End'})
    # +1 end of any regions of length 0 (Start == End)
    strling_df.loc[strling_df['Start'] == strling_df['End'], 'End'] += 1

    trf_df = parse_bed(args.trf)

    cov_pr = pr.read_bed(args.cov)

    # Annotate strling calls with corresponding PacBio calls
    strling_df = match_variants(strling_df, trf_df, args.slop)

    # Annotate strling calls with coverage overlap
    strling_df = annotate_cov(strling_df, cov_pr)

    strling_df.to_csv(args.out, sep='\t', index=False)

if __name__ == '__main__':
    main()
