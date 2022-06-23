#!/usr/bin/env python3
#Functions for manipulating fasta files in various ways

import os
import argparse
import re

def split_fasta_into_seqs(align_file, out_dir, out_file_locus_name):
    """
    Function that splits a multiple sequence alignment file in fasta format \
    into individual sequences in individual files.
    """
    with open(align_file, 'r') as raw_seqs:

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
                cluster_name = out_file_locus_name
                full_fasta_name = out_dir + '/' + cluster_name + "_" + name
                open_file = open(full_fasta_name, "w")
                open_file.write(">")
                open_file.write(name_and_seq)


def split_multiple_fastas_into_seqs(align_dir, out_dir, align_suffix):
    """
    Uses the split_fasta_into_seqs function to split up many fasta files into individual sequences. \
    Ideal for use when analyzing the bases of a concatenated alignment updated with Extensiphy.
    """

    if os.path.isdir(out_dir) == False:
        print("Specified outdir doesnt exist yet. Creating directory.")
        os.mkdir(out_dir)
    else:
        print("Specified outdir exists. Moving on.")

    path = os.path.abspath(align_dir)

    list_of_files = os.listdir(path)

    for file in list_of_files:
        path_to_align = path + '/' + file
        if os.path.isfile(path_to_align):
            if file.endswith(align_suffix):
                suffix_removed_align_name = file.strip(align_suffix)
                # print(suffix_removed_align_name)

                split_fasta_into_seqs(path_to_align, out_dir, suffix_removed_align_name)
            else:
                print("FILE NOT ENDING IN SUFFIX: ", file)

        else:
            print("DIRECTORY: ", path_to_align)
