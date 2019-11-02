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
    #parser = argparse.ArgumentParser(description = 'Produce multiple single taxon fasta files')
    parser.add_argument('-m', action='store_true', help='Turn on multi-fasta output option')
    #parser = argparse.ArgumentParser(description = 'Set function to produce a single reference fasta file')
    parser.add_argument('-s', action='store_true', help='Turn on single taxon fasta output option')
    return parser.parse_args()


def main():
    args = parse_args()

    if args.s == True:

        with open(args.align_file, 'r') as raw_seqs:

            read_seq = raw_seqs.read()
            split_seqs = read_seq.split(">")

            new_file = open(args.out_file,'w')
            new_file.write(">")
            new_file.write(split_seqs[1])

        raw_seqs.close()


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
