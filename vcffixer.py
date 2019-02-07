#!/usr/bin/env python
# Program handles an issue found in VCF files that use the reference nucleotide
# in the flattened consensus sequence when the reference is an 'N'
# instead of the nucleotide found in the reads
# USE: vcffixer.py --vcf_file [VCF_FILE] --align_file [ALIGNMENT_FILE]

import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_file')
    #parser.add_argument('--align_file')
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.vcf_file) as vcf:
        for line in vcf:
            splitter = line.split()
            ref = splitter[3:4]
            alt = splitter[4:5]
            if 'N' in ref:
                print(ref)
                print(alt)
                str_alt = ''.join(alt)
                split_alt = str_alt.split(',')
                print(split_alt[0])
                #select_alt = alt.split(',')
                #print(select_alt)

if __name__ == '__main__':
    main()
