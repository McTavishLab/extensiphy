#!/usr/bin/env python
# Program finds and fixes any off-by-one errors in an alignment

import os
import argparse
from collections import defaultdict


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

    for chunk in split_file:
        name_and_seq = ''.join(chunk)
        split_name_and_seq = name_and_seq.split('\n')
        tax_name = split_name_and_seq[0:1]
        seq = split_name_and_seq[1:2]
        # print(tax_name)
        # print(seq)
        str_name = ''.join(tax_name)
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
    print(nuc_list)
    print(len_count)
    print(len(nuc_list))
        # seq_count+=1
    #     individ_seq_count = 0
    #     seperate_seq_dict = {}
    #     sep_seq_list = (list(str_seq))
    #     for item in sep_seq_list:
    #         individ_seq_count+=1
    #         seperate_seq_dict[individ_seq_count] = item
    #     seqs_dict[str_name] = seperate_seq_dict
    #
    # for taxon, seq in seqs_dict.items():
    #
    #     # print(taxon)
    #     for key in seq:
    #         #new_list + str(key) = []
    #         nuc_list.append(seq[key])
    # print(nuc_list)
            # nuc_list.append(new_list[key])
            # print(str(key) + '----', seq[key])





if __name__ == '__main__':
    main()
