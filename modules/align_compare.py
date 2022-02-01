#!/usr/bin/env python3
# Program for comparing two sequence alignments with identical taxon labels
# Reports which positions between two identical taxon labels have non-identical bases
# OUTPUT: A .csv file recording positions and the bases contained by each alignment
# USE:
import os
import subprocess
import argparse
import pandas as pd

import sys
if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

_DEBUG = 0
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_1')
    parser.add_argument('--align_2')
    parser.add_argument('--out_file', default=None)
    return parser.parse_args()


def main():
    args = parse_args()

    # Convert each alignment into a list of lists,
    # containing split taxon names and sequences
    open_align_one = process_alignment(args.align_1)
    open_align_two = process_alignment(args.align_2)

    output = compare_and_output(open_align_one, open_align_two)

    if args.out_file != None:
        output.to_csv(args.out_file)



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

    # loop over lists of taxon names and sequences, splitting on
    # the first new line characters and take the first chunk as the taxon name
    for chunk in split_seqs:

        # Make split on first newline char, separating name and sequence
        split_name_and_seq = chunk.split('\n', 1)

        if _DEBUG:
            print(len(split_name_and_seq))

        # Filter if-statement to remove empty segments from split
        if len(split_name_and_seq) > 1:

            # add name and sequence to output dictionary
            name = split_name_and_seq[0]
            seq = split_name_and_seq[1]
            output[name] = seq


    return output
    # assert((set(seq.upper()))==set('ATGCN')), set(seq.upper())

def make_comparison(seq_1, seq_2, taxon_name):
    """Compare both sequences, noting difference.
    Return the position of the differing nucleotides and the bases themselves"""
    output = []

    # Make each seq a list and strip newline characters so theyre contiguous
    list_1 = list(seq_1.strip('\n'))
    list_2 = list(seq_2.strip('\n'))

    # Zip the lists together, base by base, position by position
    zipped_lists = list(zip(list_1, list_2))

    # Loop over the positions in the combined list, counting each position
    for num, nuc_pairs in enumerate(zipped_lists):

        # Initialize output list for a position
        position_outputs = []

        # Convert each nucleotide to uppercase
        nuc_1 = nuc_pairs[0].upper()
        nuc_2 = nuc_pairs[1].upper()

        # Check if the bases are identical
        if nuc_1 != nuc_2:

            if _DEBUG:
                print(num)
                print(nuc_pairs)

            # If not identical bases, append information to output list
            position_outputs.append(taxon_name)
            position_outputs.append(num)
            position_outputs.append(nuc_pairs[0])
            position_outputs.append(nuc_pairs[1])

        # If the output list isn't empy, added it to the list being returned
        # by the function
        if len(position_outputs) > 0:
            output.append(position_outputs)

    return output


def compare_and_output(dict_of_seqs_1, dict_of_seqs_2):
    """Loop over each entry in the first dict of seqs, identify the taxon name
    and find the matching entry in the second dict. Compare the sequences"""

    # Initialize list of lists to be converted into a df
    # Make first list the column names
    output = [['taxon', 'position', 'file_one_base', 'file_two_base']]

    # Loop over items in the dictionary from file 1
    for key, value in dict_of_seqs_1.items():
        if _DEBUG:
            print(key)

        # Find the matching taxon name and sequence in the dictionary of file 2
        opposing_seq = dict_of_seqs_2[key]
        # print(opposing_seq)

        # Pass the two sequences and the taxon name to the make_comparison function
        seq_comparison = make_comparison(value, opposing_seq, key)

        # check if the outputs of make_comparison are empty
        if len(seq_comparison) > 0:

            # If not empty, loop over contents and add the outputs
            # to the list of lists prior to df conversion
            for item in seq_comparison:
                output.append(item)

    if _DEBUG:
        print(output)

    # Convert list of lists to df, using the first list as column names
    df = pd.DataFrame(output[1:], columns=output[0])

    if _DEBUG:
        print(df)

    return df

if __name__ == '__main__':
    main()
