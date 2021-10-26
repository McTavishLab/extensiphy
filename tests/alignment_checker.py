#! /usr/bin/env python3
import argparse
import os
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_file')
    return parser.parse_args()

def main():
    args = parse_args()

    input_file = os.path.realpath(args.align_file)

    with open(input_file, 'r') as alignment:
        for line in alignment:
            print(line)

def check_alignment_lines(line):
    find_name = ">\w+"
    find_seq = "[acgt]+"

    compiled_find_name = re.compile(find_name)
    compiled_find_seq = re.compile(find_seq)

    search_name = re.search(compiled_find_name, line)
    search_seq = re.search(compiled_find_seq, line)

    




if __name__ == '__main__':
    main()
