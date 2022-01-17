#!/usr/bin/env python
# Program handles an issue found in VCF files that use the reference nucleotide
# in the flattened consensus sequence when the reference is an 'N'
# instead of the nucleotide found in the reads
# USE: vcffixer.py --vcf_file [VCF_FILE] --align_file [ALIGNMENT_FILE]
# TODO add functionality to handle dropped nucleotides in the sequence
# to prevent off-by-one errors

import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_file')
    parser.add_argument('--align_file')
    parser.add_argument('--out_file')
    return parser.parse_args()


def main():
    args = parse_args()

    nuc_set = ['A', 'C', 'G', 'T']

    # create an empty list that will hold the fixed sequence
    ref_sequence_list = []

    # initialize a count for nucleotides in the sequence
    nuc_count = 0

    # open sequence file that will have the nucleotides replaced
    seq_file = open(args.align_file, 'r')

    # read sequence file, split the file on new line characters
    # and take the first chunk as the taxon name
    read_seq = seq_file.read()
    split_seqs = read_seq.split('\n')
    tax_name = split_seqs[0:1]
    str_tax_name = ''.join(tax_name)

    # seperate out the original sequence produced by mpileup
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
    pos_duplicate_check = 0
    current_position_no_duplicates = 0
    pos_duplicate_count = 0

    with open(args.vcf_file) as vcf:
        for line_num, line in enumerate(vcf):
            if not line.startswith("#"):
                splitter = line.split()
                ref_name = splitter[0:1]
                ref_str = ''.join(ref_name)
                pos = splitter[1:2]
                if pos_duplicate_check != pos:
                    pos_duplicate_check = pos

                elif pos_duplicate_check == pos:
                    pos_duplicate_count+=1

                ref = splitter[3:4]
                alt = splitter[4:5]
                if 'N' in ref:
                    if alt in nuc_set:
                        pos_str = ''.join(pos)
                        str_alt = ''.join(alt)
                        split_alt = str_alt.split(',')
                        alt_nuc = split_alt[0]
                        replace_nuc_dict[pos_str] = alt_nuc
                    elif alt not in nuc_set:
                        pos_str = ''.join(pos)
                        str_alt = ''.join(alt)
                        split_alt = str_alt.split(',')
                        alt_nuc = 'N'
                        replace_nuc_dict[pos_str] = alt_nuc
                    else:
                        print('YOU HAVE A LARGE PROBLEM WITH YOUR VCFFIXER.PY')
                # else:
                #     ref_str = ''.join(ref)
                #     split_multiple_nucs = list(ref_str)
                #     ref_nuc_count = 0
                #     for nuc in split_multiple_nucs:
                #         ref_nuc_count+=1
                #     if ref_nuc_count > 1:
                #         if pos_duplicate_count == 1:
                #             pos_str = ''.join(pos)
                #             split_alt = ref_str.split(',')
                #             alt_nuc = split_alt[0]
                #             replace_nuc_dict[pos_str] = alt_nuc
                #
                #         elif pos_duplicate_count > 1:
                #             pos_str = ''.join(pos)
                #             str_alt = ''.join(alt)
                #             split_alt = str_alt.split(',')
                #             alt_nuc = split_alt[0]
                #             replace_nuc_dict[pos_str] = alt_nuc





    # loop through the dictionary and find the positions in the list
    # that matches the position for the new nucleotide in the dictionary
    # replace the nucleotide in the list at that position
    # with the nucleotide from the dictionary
    for key, value in replace_nuc_dict.items():
        print(ref_sequence_list[int(key) - 1])
        print(key)
        print(value)
        # print("~~~~~~~~~~~~~~~~~~~~~~~")

    #########################################################
    #original code
        ref_sequence_list[int(key) - 1] = value
        # print(ref_sequence_list[int(key)])



    #########################################################
    #Original output file code
    #
    # # turn the list into a string
    # fixed_seq = ''.join(ref_sequence_list)
    #
    # # write the taxon name, a new line character and the fixed sequence to file
    # new_file = open(args.out_file,'w')
    # new_file.write(str_tax_name)
    # new_file.write('\n')
    # new_file.write(fixed_seq)

if __name__ == '__main__':
    main()
