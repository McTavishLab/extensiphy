#!/usr/bin/env python
# Program finds and fixes any off-by-one errors in an alignment

import os
import argparse
from collections import defaultdict
import itertools


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_file')
    # parser.add_argument('--out_file')
    return parser.parse_args()

def main():
    args = parse_args()

    # with open(args.align_file) as seq_file:
    #     for line in seq_file:
    seqs_dict = {}
    nuc_list = []
    seq_file = open(args.align_file, 'r')
    read_seq = seq_file.read()
    split_file = read_seq.split('>')
    seq_count = 0
    str_split_file = ''.join(split_file[1])
    split_file_sequence = str_split_file.split('\n')
    # print(split_file_sequence[1])
    for letter in split_file_sequence[1]:
        nuc_list.append([])
    # print(nuc_list)

    name_list = {}
    for position, chunk in enumerate(split_file):
        name_and_seq = ''.join(chunk)
        split_name_and_seq = name_and_seq.split('\n')
        tax_name = split_name_and_seq[0:1]
        seq = split_name_and_seq[1:2]
        # print(tax_name)
        # print(seq)
        str_name = ''.join(tax_name)
        name_list[position] = str_name
        str_seq = ''.join(seq)
        # print(str_seq)
        len_count = 0
        # for letter in str_seq:
        #     nuc_list.append([])
        for num1, letter in enumerate(str_seq):
            for num_list, list in enumerate(nuc_list):
                if num1 == num_list:
                    len_count+=1
                    list.append(letter)
    # print(nuc_list)
    # print(len_count)
    # print(len(nuc_list))

    list_of_identical_nucleotides = 0
    list_of_non_identical_nucleotides = 0
    position_of_SNPs_start = ''
    taxon_snp_dict = {}
    first_list = nuc_list[0]

    # create the dictionary of taxa that are indicated by their position in the alignment
    # This will be updated whenever a nucleotide at the correct position is a snp
    for position in range(len(first_list)):
        taxon_snp_dict[position] = 0

    # loop through the list of lists containing each nucleotide at the same position
    # for all taxa
    # start the counters for identical sites and SNPs at the beginning of each loop
    for position, list in enumerate(nuc_list):
        identical_nucleotides = 0
        non_identical_nucleotides = 0
        first_nucleotide = ''

        # loop through each position in the SNP position dictionary
        for position in taxon_snp_dict:

            # Loop through each letter in the list of nucleotides
            # establishing the first letter as the comparison letter for the rest of the identical_nucleotides
            # TODO implement majority rule for the first, consensus nucleotide
            for taxon, letter in enumerate(list):
                letter = letter.upper()
                if len(first_nucleotide) == 0:
                    first_nucleotide = letter
                elif len(first_nucleotide) != 0:
                    other_letter = letter

                    # Establish if the next letters are identical to the first letter
                    # and update the appropriate counter
                    if other_letter == first_nucleotide:
                        identical_nucleotides+=1
                    elif other_letter != first_nucleotide:
                        non_identical_nucleotides+=1
                        if position == taxon:
                            taxon_snp_dict[position]+=1

    final_dict = {}
    for position, name in name_list.items():
        for pos, count in taxon_snp_dict.items():
            if position == pos:
                final_dict[name] = count
    print(final_dict)






if __name__ == '__main__':
    main()
