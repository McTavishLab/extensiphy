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
    for chunk in split_file:
        name_and_seq = ''.join(chunk)
        split_name_and_seq = name_and_seq.split('\n')
        #print(split_name_and_seq[1:2])
        tax_name = split_name_and_seq[0:1]
        seq = split_name_and_seq[1:2]
        # print(tax_name)
        # print(seq)
        str_name = ''.join(tax_name)
        str_seq = ''.join(seq)
        seq_count+=1
        individ_seq_count = 0
        seperate_seq_dict = {}
        sep_seq_list = (list(str_seq))
        for item in sep_seq_list:
            individ_seq_count+=1
            seperate_seq_dict[individ_seq_count] = item
        seqs_dict[str_name] = seperate_seq_dict
        # nuc_list.append(seperate_seq_dict.copy())
    print(seqs_dict)
        # seqs_dict[str_name] = str_seq
    # print(seqs_dict)
    # taxon_count = 0
    # for key, value in seqs_dict.items():
    #     split_list = value.split()
    #     print(list(split_list))
    # for key, value in zip(seqs_dict.items()):
    #     print(value, value)
        # for item in split_name_and_seq:
        #     print(item[0:4])
            # seqs_dict[item[0:1]] = item[1:2]




if __name__ == '__main__':
    main()
