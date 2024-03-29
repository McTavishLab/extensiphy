#!/usr/bin/env python3

import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description = 'Replace degenerate nucleotides in a fasta file with Ns')
    parser.add_argument('--align_file')
    parser.add_argument('--output')
    return parser.parse_args()


def main():
    args = parse_args()

    degen_nucs = {'-'}

    output_data = []

    data = args.align_file

    read_data = open(data,'r').read()

    split_data = read_data.split('>')

    # print(split_data)
    for num, chunk in enumerate(split_data):
        if len(chunk) == 0:
            del split_data[num]


    for chunk in split_data:
        split_newlines = chunk.split('\n', 1)
        # print(split_newlines)

        # Make sure there isn't a list of 1 empty string
        if len(split_newlines) > 1:
            # print("SPLIT FOUND")
            fixed_seq = broken_seq_handler(split_newlines)
            # print(fixed_seq)

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


def broken_seq_handler(seq_list):
    # print(seq_list)

    returned_seq = []
    #Make sure theres actually sequence information in the seq_list variable
    if len(seq_list) >= 2:
        # print(seq_list)
        seq_id = seq_list[0]
        joined_seq = ''.join(seq_list[1:])

        gapless_joined_seq = remove_gaps(joined_seq)

        gapless_joined_seq = gapless_joined_seq.replace('\n','')

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



if __name__ == '__main__':
    main()
