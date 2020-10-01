#! /usr/bin/env python
# program produces the best reference file from the original alignment
# currently, this program just isolates the first sequence in the file
# this handles if there are new line breaks and the sequence isn't read as a single line (some genome files)
import os
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_file')
    parser.add_argument('--out_file')
    parser.add_argument('--ref_select')
    #parser = argparse.ArgumentParser(description = 'Produce multiple single taxon fasta files')
    parser.add_argument('-m', action='store_true', help='Turn on multi-fasta output option')
    #parser = argparse.ArgumentParser(description = 'Set function to produce a single reference fasta file')
    parser.add_argument('-r', action='store_true', help='Turn on single taxon RANDOM reference fasta output option')
    parser.add_argument('-s', action='store_true', help='Turn on single taxon SEARCH reference fasta output option')
    return parser.parse_args()


def main():
    args = parse_args()

    if args.r == True:

        with open(args.align_file, 'r') as raw_seqs:

            read_seq = raw_seqs.read()
            split_seqs = read_seq.split(">")

            new_file = open(args.out_file,'w')
            new_file.write(">")
            new_file.write(split_seqs[1])

        raw_seqs.close()

    elif args.s == True:

        with open(args.align_file, 'r') as raw_seqs:

            read_seq = raw_seqs.read()
            split_seqs = read_seq.split(">")
            # print(split_seqs)
            # Construct regex of name of reference selected by the user
            ref_name = args.ref_select
            # ref_name = ref_name + "\n"
            ref_name = "^" + ref_name + "\n"
            # print(ref_name)
            ref_name_compile = re.compile(ref_name)
            
            # Loop through sequences and identify which one has the desired reference
            for taxon in split_seqs:
                ref_search = re.search(ref_name_compile, taxon)
                if ref_search is not None:
                    #print(taxon)
                    ref_seq_file = open(args.out_file,'w')
                    ref_seq_file.write(">")
                    ref_seq_file.write(taxon)
                    ref_seq_file.close()
                #elif ref_search is None:
                #    try:
                #        raise ValueError('DESIRED REFERENCE NOT FOUND IN ALIGNMENT FILE, CHECK YOUR SPELLING AND MAKE SURE THE REFERENCE IS IN THE FILE')
                        #print("DESIRED REFERENCE NOT FOUND IN ALIGNMENT FILE, CHECK YOUR SPELLING AND MAKE SURE THE REFERENCE IS IN THE FILE")

    elif args.m == True:

        with open(args.align_file, 'r') as raw_seqs:

            read_seq = raw_seqs.read()
            split_seqs = read_seq.split(">")
            #print(split_seqs[1])
            for name_and_seq in split_seqs:
                if len(name_and_seq) > 2: 
                    split_name = name_and_seq.split("\n", 1)
                    #print(split_name)
                    name = split_name[0]
                    #print(name)
                    seq = split_name[1]
                    cluster_name = args.out_file
                    full_fasta_name = cluster_name + "_" + name
                    open_file = open(full_fasta_name, "w")
                    open_file.write(">")
                    open_file.write(name_and_seq)
            
            #open_file.write(">")
            #open_file.write(name)
            #open_file.write(seq)
            #open_file.close()


            #for item in split_seqs[1]:
            #    split_name_and_seq = item.split("\n")
            #    print(split_name_and_seq)

            
            #new_file = open(args.out_file,'w')
            #new_file.write(">")
            #new_file.write(split_seqs[1])


if __name__ == '__main__':
    main()
