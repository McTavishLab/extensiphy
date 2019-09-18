#! /usr/bin/env python
# program produces the best reference file from the original alignment
# currently, this program just isolates the first sequence in the file
# this handles if there are new line breaks and the sequence isn't read as a single line (some genome files)
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_file')
    parser.add_argument('--out_file')
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.align_file, 'r') as raw_seqs:

        read_seq = raw_seqs.read()
        split_seqs = read_seq.split(">")

        new_file = open(args.out_file,'w')
        new_file.write(">")
        new_file.write(split_seqs[1])

    raw_seqs.close()


if __name__ == '__main__':
    main()
