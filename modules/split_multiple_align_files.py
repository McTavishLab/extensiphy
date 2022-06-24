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
    parser.add_argument('--genome_dir', help='Directory of complete genomes')
    # parser.add_argument('--ref_select')
    # parser.add_argument('-m', action='store_true', help='Turn on multi-fasta output option')
    # parser.add_argument('-r', action='store_true', help='Turn on single taxon RANDOM reference fasta output option')
    # parser.add_argument('-s', action='store_true', help='Turn on single taxon SEARCH reference fasta output option')
    return parser.parse_args()


def main():
    args = parse_args()

    if os.path.isdir(args.out_dir) == False:
        print("Specified outdir doesnt exist yet. Creating directory.")
        os.mkdir(args.out_dir)
    else:
        print("Specified outdir exists. Moving on.")

    # if os.path.isdir(args.out_dir + '/contiguous_genome_files') == False:
    #     os.chdir(args.out_dir)
    #     # os.mkdir('contiguous_genome_files')
    #     os.mkdir('split_locus_files')
    #
    # else:
    os.chdir(args.out_dir)

    working_dir = os.path.abspath(os.getcwd())

    matching_seqs_dir = working_dir + '/match_outputs'
    os.mkdir(matching_seqs_dir)

    make_genomes_contiguous(working_dir, args.genome_dir, args.align_suffix)

    split_multiple_fastas_into_seqs(working_dir, args.align_dir, args.align_suffix)

    match_long_with_loci(working_dir + '/split_locus_files', working_dir + '/contiguous_genome_files', matching_seqs_dir)




if __name__ == '__main__':
    main()