#!/usr/bin/env python

import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description = 'Replace degenerate nucleotides in a fasta file with Ns')
    parser.add_argument('--align_file')
    parser.add_argument('--output', default='ref_nogaps.fas')
    return parser.parse_args()


def broken_seq_handler(seq_list):

    returned_seq = []

    # print(seq_list)

    # print("waffle1")
    #Make sure theres actually sequence information in the seq_list variable
    assert len(seq_list[0]) > 0 and len(seq_list[1]) > 0

    # if len(seq_list) > 2:
    #     # print("waffle2")
    #     # print(seq_list)
    seq_id = seq_list[0]
    joined_seq = ''.join(seq_list[1:])

    gapless_joined_seq = remove_gaps(joined_seq)

    #Test to make sure we removed all new lines from the sequence
    assert '\n' not in gapless_joined_seq

    # Test to make sure we've removed gap characters
    assert '-' not in gapless_joined_seq

    returned_seq.append(seq_id)
    returned_seq.append(gapless_joined_seq)

    assert len(returned_seq) == 2

    return returned_seq


def remove_gaps(seq):
    fixed_seq = seq.replace('-','')
    return fixed_seq

def main():
    args = parse_args()

    degen_nucs = {'-'}

    output_data = []

    data = args.align_file

    read_data = open(data,'r').read()

    split_data = read_data.split('>')

    for chunk in split_data:
        split_newlines = chunk.split('\n')

        # Make sure there isn't a list of 1 empty string
        if len(split_newlines) > 1:
            fixed_seq = broken_seq_handler(split_newlines)

            output_data.append(fixed_seq)

    output_file = open(args.output, 'w')
    # print(len(output_data))
    for name_and_seq in output_data:
        if name_and_seq != None:
            seq_id = name_and_seq[0]
            seq = name_and_seq[1]
            # print(seq_id)
            # print(seq)

            output_file.write('>')
            output_file.write(seq_id)
            output_file.write('\n')
            output_file.write(seq)
            output_file.write('\n')

    output_file.close()




if __name__ == '__main__':
    main()
