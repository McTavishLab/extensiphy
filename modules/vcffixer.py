#!/usr/bin/env python
# Program handles an issue found in VCF files that use the reference nucleotide
# in the flattened consensus sequence when the reference is an 'N'
# instead of the nucleotide found in the reads
# USE: vcffixer.py --vcf_file [VCF_FILE] --align_file [ALIGNMENT_FILE]
# TODO add functionality to handle dropped nucleotides in the sequence
# to prevent off-by-one errors

import os
import subprocess
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
    ref_sequence_dict = {}

    # initialize dictionary for location of nucleotide to be replaced
    # and the nucleotide that will go there
    replace_nuc_dict = {}

    # initialize a count for nucleotides in the sequence
    nuc_count = 0

    # begin looping through VCF file to look for positions that will be replaced
    # and add them to the dictionary
    pos_duplicate_check = 0
    current_position_no_duplicates = 0
    pos_duplicate_count = 0

    # initialize reference name variable. this may not be necessary later
    ref_str = ''

    proc_seq = process_alignment(args.align_file, ref_sequence_list, nuc_count, ref_str)

    last_line = subprocess.check_output(['tail', '-1', args.vcf_file]).decode('UTF-8')
    split_last_line = last_line.split()
    vcf_length = split_last_line[1]
    print(vcf_length)

    process_vcf(args.vcf_file, replace_nuc_dict)

    # # open sequence file that will have the nucleotides replaced
    # seq_file = open(args.align_file, 'r')
    #
    # # read sequence file, split the file on new line characters
    # # and take the first chunk as the taxon name
    # read_seq = seq_file.read()
    # split_seqs = read_seq.split('\n')
    # # tax_name = split_seqs[0:1]
    # tax_name = split_seqs[0]
    # print(tax_name)
    # str_tax_name = ''.join(tax_name)
    #
    # # seperate out the original sequence produced by mpileup
    # seq = split_seqs[1:2]
    # seq_str = seq[0]
    # for nuc in seq_str:
    #     ref_sequence_list.append(nuc)
    #     nuc_count+=1


    # # initialize reference name variable. this may not be necessary later
    # ref_str = ''

    # # initialize dictionary for location of nucleotide to be replaced
    # # and the nucleotide that will go there
    # replace_nuc_dict = {}
    #
    # # begin looping through VCF file to look for positions that will be replaced
    # # and add them to the dictionary
    # pos_duplicate_check = 0
    # current_position_no_duplicates = 0
    # pos_duplicate_count = 0

    # last_line = subprocess.check_output(['tail', '-1', args.vcf_file]).decode('UTF-8')
    # split_last_line = last_line.split()
    # vcf_length = split_last_line[1]
    # print(vcf_length)


    # with open(args.vcf_file) as vcf:
    #     for line_num, line in enumerate(vcf):
    #         if not line.startswith("#"):
    #             splitter = line.split()
    #             # ref_name = splitter[0:1]
    #             ref_name = splitter[0]
    #             # print(ref_name)
    #             ref_str = ''.join(ref_name)
    #             # pos = splitter[1:2]
    #             pos = splitter[1]
    #             # print(pos)
                # if pos_duplicate_check != pos:
                #     pos_duplicate_check = pos
                #
                # elif pos_duplicate_check == pos:
                #     pos_duplicate_count+=1
                #
                # ref = splitter[3:4]
                # alt = splitter[4:5]
                # if 'N' in ref:
                #     if alt in nuc_set:
                #         pos_str = ''.join(pos)
                #         str_alt = ''.join(alt)
                #         split_alt = str_alt.split(',')
                #         alt_nuc = split_alt[0]
                #         replace_nuc_dict[pos_str] = alt_nuc
                #     elif alt not in nuc_set:
                #         pos_str = ''.join(pos)
                #         str_alt = ''.join(alt)
                #         split_alt = str_alt.split(',')
                #         alt_nuc = 'N'
                #         replace_nuc_dict[pos_str] = alt_nuc
                #     else:
                #         print('YOU HAVE A LARGE PROBLEM WITH YOUR VCFFIXER.PY')





    # loop through the dictionary and find the positions in the list
    # that matches the position for the new nucleotide in the dictionary
    # replace the nucleotide in the list at that position
    # with the nucleotide from the dictionary
    # print(len(replace_nuc_dict))
    # print(replace_nuc_dict)

    # print(replace_nuc_dict[vcf_length])
    # find_len_diff = find_length_difference(replace_nuc_dict, vcf_length)
    # print(find_len_diff)
    # find_first_position(replace_nuc_dict, vcf_length)
    # for key, value in replace_nuc_dict.items():
    #     ref_sequence_list[int(key) - 1] = value
    #
    # # turn the list into a string
    # fixed_seq = ''.join(ref_sequence_list)
    #
    # # write the taxon name, a new line character and the fixed sequence to file
    # new_file = open(args.out_file,'w')
    # new_file.write(str_tax_name)
    # new_file.write('\n')
    # new_file.write(fixed_seq)

def process_vcf(vcf_file, dict):
    """Read VCF file and add nucleotides and positions to a dictionary. \
    Dict protects against multiple nucleotide options at a single position"""
    with open(vcf_file) as vcf:
        for line_num, line in enumerate(vcf):
            if not line.startswith("#"):
                splitter = line.split()
                ref_name = splitter[0]
                ref_str = ''.join(ref_name)
                pos = splitter[1]
                # ref = splitter[3:4]
                # alt = splitter[4:5]
                ref = splitter[3]
                alt = splitter[4]
                # print(ref)
                # print(alt)

                dupe_check = check_dict_duplicate(dict, pos)
                if dupe_check == False:
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



def check_dict_duplicate(dict, check_key):
    """Check if a key/value is already present in a dictionary."""
    if check_key in dict.keys():
        print("DUPLICATE POSITION FOUND at position: ", check_key)
        return False

def process_alignment(align_file, output_dict, base_count, taxon_name):
    """read alignment file, add individual nucleotides to a list, count the bases and store the taxon name"""
    print("Processing alignment.")
    # open sequence file that will have the nucleotides replaced
    seq_file = open(align_file, 'r')

    # read sequence file, split the file on new line characters
    # and take the first chunk as the taxon name
    read_seq = seq_file.read()
    split_seqs = read_seq.split('\n')
    # tax_name = split_seqs[0:1]
    tax_name = split_seqs[0]
    # print(tax_name)
    str_tax_name = ''.join(tax_name)

    # seperate out the original sequence produced by mpileup
    seq = split_seqs[1]
    # seq_str = seq[0]
    # print(seq_str)
    for nuc in seq:
        output_dict.append(nuc)
        base_count+=1

    print("Alignment nucleotide count: ", base_count)


def find_length_difference(dict_of_positions, vcf_length):
    int_len = int(vcf_length)
    str_vcf_length = str(int_len)
    if str_vcf_length not in dict_of_positions.keys():
        find_length_difference(dict_of_positions, (int_len - 1))
    elif str_vcf_length in dict_of_positions.keys():
        return str_vcf_length

def find_first_position(dict_of_positions, end_len):
    starting_pos = 0
    while starting_pos < int(end_len) and starting_pos < 800:
        if str(starting_pos) in dict_of_positions.keys():
            print("--")
        else:
            starting_pos+=1
            print("NOT FOUND: ", starting_pos)


if __name__ == '__main__':
    main()
