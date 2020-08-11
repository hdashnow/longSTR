# Parse TRF ngs dat file into tab delimited file
# Assumes TRF was run on aligned contigs.
# Each @ is a contig labeled with it's mapping position,
# each line is an STR

import argparse
import sys
import os
from strtools import normalise_str
import pandas as pd
#import numpy as np
import pyranges as pr

def parse_args():
    # top-level parser
    parser = argparse.ArgumentParser()

    parser.add_argument('--strling', type=str,
                        help='STRling variants')
    parser.add_argument('--trf', type=str, nargs='+',
                        help='STR variants from TRF in bed format (two, one for each haplotype)')
    parser.add_argument('--cov', type=str, nargs='+',
                        help='Contig coverage in bed format (two, one for each haplotype)')
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

def match_variants(strling_df, trf_hap1_df, trf_hap2_df, slop=100):
    """Match up pacbio trf variants to their corresponding strling variants
        - Within X bp of slop
        - Same repeat unit
        Any leftover pacbio trf variants added on to the end"""
    strling_pr = pr.PyRanges(strling_df)
    trf_pr = pr.PyRanges(trf_hap1_df)

    # Note, this will annotate with the closest physical locus
    # Could break down by repeat unit first then do this?
    nearest_pr = strling_pr.nearest(trf_pr)
    #nearest_pr.to_csv('tmp.tsv', sep='\t')
    nearest_df = nearest_pr.df

    # Remove pacbio variants more than slop bp away
    nearest_columns = ['Start_b', 'End_b', 'repeatunit_b', 'period', 'length_ru', 'length_bp', 'indel', 'sample_b', 'Distance']
    nearest_df.loc[nearest_df.Distance > slop, nearest_columns] = None #np.nan
    #nearest_df.to_csv('tmp2.tsv', sep='\t')

    # Remove pacbio variants with a different repeat unit
    diff_ru = nearest_df['repeatunit'].apply(normalise_str) != nearest_df['repeatunit_b'].apply(normalise_str)
    nearest_df.loc[diff_ru, nearest_columns] = None #np.nan
    nearest_df.to_csv('tmp2.tsv', sep='\t')

def annotate_cov(strling_df, cov_hap1_pr, cov_hap2_pr):
    """Annotate strling calls with contig coverage overlap"""
    strling_pr = pr.PyRanges(strling_df)

    for i, hap_pr in zip((1,2), (cov_hap1_pr, cov_hap2_pr)):
        cov_pr = strling_pr.coverage(hap_pr)
        cov_df = cov_pr.df
        # If FractionOverlaps < 1, set NumberOverlaps to 0 
        cov_df['NumberOverlaps'][cov_df.FractionOverlaps < 1] = 0
        strling_df[f'Hap{i}Cov'] = cov_df['NumberOverlaps']

    return(strling_df)

def main():
    args = parse_args()
    if len(args.trf) != 2:
        exit(f'ERROR: Expected 2 trf files, got {len(args.trf)}: {args.trf}')
    if len(args.cov) != 2:
        exit(f'ERROR: Expected 2 cov files, got {len(args.cov)}: {args.cov}')

    # Parse inputs as pandas data frame (df) or pyranges (pr) objects
    strling_df = parse_bed(args.strling)

    trf_hap1_df = parse_bed(args.trf[0])
#    trf_hap1_pr = pr.PyRanges(trf_hap1_df)

    trf_hap2_df = parse_bed(args.trf[1])
#    trf_hap2_pr = pr.PyRanges(trf_hap2_df)

    cov_hap1_pr = pr.read_bed(args.cov[0])
    cov_hap2_pr = pr.read_bed(args.cov[1])

    # Annotate strling calls with corresponding PacBio calls
    match_variants(strling_df, trf_hap1_df, trf_hap2_df)
    # Annotate strling calls with coverage overlap
    #strling_df = annotate_cov(strling_df, cov_hap1_pr, cov_hap2_pr)

    #print(strling_df)

if __name__ == '__main__':
    main()
