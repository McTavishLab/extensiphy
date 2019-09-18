#! /usr/bin/python

import os
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_file')
    parser.add_argument('--out_file')
    parser.add_argument('--locus_id')
    return parser.parse_args()


def main():
    args = parse_args()
    
    taxa_count = 0
    use_lines = False
    full_sequences = {}
    current_seq = ''
    stored_name = ''
    locus_id = "(>1:\d+-\d+\s+(\+|\-)\s" + args.locus_id + "\s(.|\n)+?)="
    locus_id_grabber = re.compile(locus_id, re.S)
    #print(locus_id) 
    alignment = open(args.align_file, "r")
    read_align = alignment.read()
    locus_search = re.search(locus_id_grabber, read_align)
    if locus_search:
        print(locus_search.group())

    




if __name__ == '__main__':
    main()
