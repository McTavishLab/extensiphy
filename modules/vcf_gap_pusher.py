#!/usr/bin/env python3

import os
import argparse
import pathlib
import shutil
import pandas as pd
from itertools import islice

def parse_args():
    parser = argparse.ArgumentParser(prog='vcf_gap_pusher.py', \
        description='move extensiphy vcf positions to the correct place.')
    parser.add_argument('--gap_csv')
    parser.add_argument('--vcf')
    parser.add_argument('--bases_csv')
    return parser.parse_args()

def main():
    args = parse_args()


    gap_df = pd.read_csv(args.gap_csv)

    # read_vcf = pd.read_table(args.vcf, sep='\s+', skiprows=23, escapechar='#')
    vcf = smart_read_vcf(args.vcf)

    interest_bases = pd.read_csv(args.bases_csv)

    print(interest_bases.columns)

    real_position_df = get_real_position(vcf, gap_df, interest_bases)



def get_real_position(vcf, gaps_df, bai_df):

    vcf_columns = vcf.columns
    print(vcf_columns)

    for idx, row in bai_df.iterrows():
        bai_pos = row['positions']
        # print(bai_pos)

        earlier_gap_rows = gaps_df.loc[gaps_df['gap_position'] <= bai_pos]
        # print(earlier_gap_rows)
        num_gaps = len(earlier_gap_rows)

        updated_vcf_position = bai_pos - num_gaps

        print('position of base of interest: ', bai_pos)
        print('position in vcf file to find that base: ', updated_vcf_position)
        print("#################################")



def smart_read_vcf(unread_vcf_obj):
    number_double_hash = 0
    with open(unread_vcf_obj) as vcf_file:
        head = list(islice(vcf_file, 50))
    for line in head:
        if "##" in line:
            number_double_hash+=1

    read_vcf = pd.read_table(unread_vcf_obj, sep='\s+', skiprows=number_double_hash, escapechar='#')

    return read_vcf



if __name__ == '__main__':
    main()
