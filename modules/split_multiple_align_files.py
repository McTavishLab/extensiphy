#!/usr/bin/env python3

import os
import argparse
import re
from fasta_manipulation import *

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_dir', help='Directory of alignments you wish to split')
    parser.add_argument('--out_dir', help='Directory you wish to put the split sequence files. \
    Make the dir if it doesnt exist already.')
    parser.add_argument('--align_suffix', default='.fasta', help='Suffix of alignment files you wish to operate on.')
    # parser.add_argument('--ref_select')
    # parser.add_argument('-m', action='store_true', help='Turn on multi-fasta output option')
    # parser.add_argument('-r', action='store_true', help='Turn on single taxon RANDOM reference fasta output option')
    # parser.add_argument('-s', action='store_true', help='Turn on single taxon SEARCH reference fasta output option')
    return parser.parse_args()


def main():
    args = parse_args()

    split_multiple_fastas_into_seqs(args.align_dir, args.out_dir, args.align_suffix)

    # if os.path.isdir(args.out_dir) == False:
    #     print("Specified outdir doesnt exist yet. Creating directory.")
    #     os.mkdir(args.out_dir)
    # else:
    #     print("Specified outdir exists. Moving on.")
    #
    # path = os.path.abspath(args.align_dir)
    #
    # list_of_files = os.listdir(path)
    #
    # for file in list_of_files:
    #     path_to_align = path + '/' + file
    #     if os.path.isfile(path_to_align):
    #         if file.endswith(args.align_suffix):
    #             suffix_removed_align_name = file.strip(args.align_suffix)
    #             print(suffix_removed_align_name)
    #
    #
    #
    #             split_fasta_into_seqs(path_to_align, args.out_dir, suffix_removed_align_name)
    #
    #     else:
    #         print("DIRECTORY: ", path_to_align)





if __name__ == '__main__':
    main()
