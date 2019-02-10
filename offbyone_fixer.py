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
    # print(nuc_list)
    # print(len_count)
    # print(len(nuc_list))

    list_of_identical_nucleotides = 0
    list_of_non_identical_nucleotides = 0
    position_of_SNPs_start = ''
    for position, list in enumerate(nuc_list):
        identical_nucleotides = 0
        non_identical_nucleotides = 0
        first_nucleotide = ''
        for letter in list:
            letter = letter.upper()
            if len(first_nucleotide) == 0:
                first_nucleotide = letter
            elif len(first_nucleotide) != 0:
                other_letter = letter
                if other_letter == first_nucleotide:
                    identical_nucleotides+=1
                elif other_letter != first_nucleotide:
                    non_identical_nucleotides+=1

        if identical_nucleotides == len(list) - 1:
            print("identical list of letters")
            print(list)
            list_of_identical_nucleotides+=1
            if list_of_identical_nucleotides == 20:
                list_of_non_identical_nucleotides = 0
                list_of_identical_nucleotides = 0

        elif non_identical_nucleotides > 1:
            print("non identical list of letters")
            print(list)
            list_of_non_identical_nucleotides+=1
            if list_of_non_identical_nucleotides == 150:
                print("pretty big problems")
                list_of_non_identical_nucleotides = 0

        # elif list_of_identical_nucleotides >= 10:
        #     list_of_non_identical_nucleotides = 0
        #     print("restarted list")
        #
        # elif list_of_non_identical_nucleotides >= 500:
        #     print("major friggin problems")




    # for list in nuc_list:
    #     nuc_chunk+=1
    #     #for i in range(len(list)):
    #     for i in list:
    #         print(i)
    # #       for j in range(i + 1, len(list)):
    #         for j in range(i + 1, len(list)):
    #             print(j)
    #             if list[i].upper() != list[j].upper() and \
    #              list[i].upper() != 'N' and list[j].upper() != 'N':
    #                 SNPs+=1
    #                 # homologous = 0
    #                 if SNPs == 100:
    #                     snp_chunk+=1
    #                     SNPs = 0
    #                     homologous = 0
    #             elif list[i].upper() == list[j].upper() and \
    #              list[i].upper() != 'N' and list[j].upper() != 'N':
    #                 homologous+=1
    #                 if homologous == 100:
    #                     hom_chunk = 1
    #                     SNPs = 0
    #                     homologous = 0
    # print(nuc_chunk)
    # print("Number of 100 nucleotide segments with more SNPs than homologous sites:")
    # #print(SNPs)
    # print(snp_chunk)
    # print("Number of 100 nucleotide segments with mor homologous sites than SNPs:")
    # #print(homologous)
    # print(hom_chunk)
    # assert(snp_chunk <= hom_chunk)






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
