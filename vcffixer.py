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
    parser.add_argument('--align_file')
    return parser.parse_args()


def main():
    args = parse_args()

    ref_sequence_list = []
    nuc_count = 0
    # open sequence file that will have the nucleotides replaced
    seq_file = open(args.align_file, 'r')
    read_seq = seq_file.read()
    split_seqs = read_seq.split('\n')
    seq = split_seqs[1:2]
    seq_str = seq[0]
    for nuc in seq_str:
        ref_sequence_list.append(nuc)
        nuc_count+=1


    # initialize reference name variable. this may not be necessary later
    ref_str = ''

    # initialize dictionary for location of nucleotide to be replaced
    # and the nucleotide that will go there
    replace_nuc_dict = {}

    # begin looping through VCF file to look for positions that will be replaced
    # and add them to the dictionary
    with open(args.vcf_file) as vcf:
        for line in vcf:
            splitter = line.split()
            ref_name = splitter[0:1]
            ref_str = ''.join(ref_name)
            pos = splitter[1:2]
            ref = splitter[3:4]
            alt = splitter[4:5]
            if 'N' in ref:
                pos_str = ''.join(pos)
                str_alt = ''.join(alt)
                split_alt = str_alt.split(',')
                alt_nuc = split_alt[0]
                replace_nuc_dict[pos_str] = alt_nuc

    for key, value in replace_nuc_dict.items():
        ref_sequence_list[int(key) - 1] = value

    fixed_seq = ''.join(ref_sequence_list)

    new_file = open('cns_fixed.fa','w')
    new_file.write('>' + ref_str)
    new_file.write('\n')
    new_file.write(fixed_seq)

if __name__ == '__main__':
    main()
