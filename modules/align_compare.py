#!/usr/bin/env python3
# Program for comparing two sequence alignments with identical taxon labels
# Reports which positions between two identical taxon labels have non-identical bases
# OUTPUT: A .csv file recording positions and the bases contained by each alignment
# USE:
import os
import subprocess
import argparse
import pandas

import sys
if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

_DEBUG = 1
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_1')
    parser.add_argument('--align_2')
    parser.add_argument('--out_file')
    return parser.parse_args()


def main():
    args = parse_args()

    # Convert each alignment into a list of lists,
    # containing split taxon names and sequences
    open_align_one = process_alignment(args.align_1)
    open_align_two = process_alignment(args.align_2)

    output = compare_and_output(open_align_one, open_align_two)





def process_alignment(align_file):
    """read consensus sequence file, add names and sequences to a dictionary"""
    if _DEBUG:
        print("Processing alignment.")
    #make list to hold split taxon names and sequences
    output = {}

    # open sequence file that will have the nucleotides replaced
    seq_file = open(align_file, 'r')

    # read sequence file, split the file on ">" characters
    read_seq = seq_file.read()
    split_seqs = read_seq.split('>')
    # print(split_seqs)
    if _DEBUG:
        print(len(split_seqs))

    # loop over lists of taxon names and sequences, splitting on new line characters
    # and take the first chunk as the taxon name
    for chunk in split_seqs:
        split_name_and_seq = chunk.split('\n', 1)
        if _DEBUG:
            print(len(split_name_and_seq))
        if len(split_name_and_seq) > 1:
            name = split_name_and_seq[0]
            seq = split_name_and_seq[1]
            output[name] = seq


    return output
    # assert((set(seq.upper()))==set('ATGCN')), set(seq.upper())

def make_comparison(seq_1, seq_2):
    """Compare both sequences, noting difference.
    Return the position of the differing nucleotides and the bases themselves"""
    output = []

    # Make each seq a list
    list_1 = list(seq_1.strip('\n'))
    list_2 = list(seq_2.strip('\n'))

    zipped_lists = list(zip(list_1, list_2))

    for num, nuc_pairs in enumerate(zipped_lists):
        position_outputs = []
        nuc_1 = nuc_pairs[0].upper()
        nuc_2 = nuc_pairs[1].upper()
        if nuc_1 != nuc_2:
            if _DEBUG:
                print(num)
                print(nuc_pairs)
            position_outputs.append(num)
            position_outputs.append(nuc_pairs[0])
            position_outputs.append(nuc_pairs[1])
        if len(position_outputs) > 0:
            output.append(position_outputs)

    return output

def format_comparisons():
    """organizes comparson data into dataframe"""


def compare_and_output(dict_of_seqs_1, dict_of_seqs_2):
    """Loop over each entry in the first dict of seqs, identify the taxon name
    and find the matching entry in the second dict. Compare the sequences"""
    output = {}

    for key, value in dict_of_seqs_1.items():
        if _DEBUG:
            print(key)
        opposing_seq = dict_of_seqs_2[key]
        # print(opposing_seq)
        seq_comparison = make_comparison(value, opposing_seq)
        if len(seq_comparison) > 0:
            output[key] = seq_comparison

    print(output)

if __name__ == '__main__':
    main()
