#! /usr/bin/python

import os
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_dir')
    parser.add_argument('--out_file')
    parser.add_argument('--len_filter')
    return parser.parse_args()

def main():
    args = parse_args()

    dir_of_aligns = args.align_dir
    
    # check if path ends with a "/"
    if dir_of_aligns.endswith("/"):
        print("Pathing acceptable")
        
    else:
        print("Fixing pathing")
        dir_of_aligns = dir_of_aligns + "/"
        #print(dir_of_aligns)

    dict_of_names_and_seqs = {}
    dict_of_loci_and_len_and_order = {}
    
    locus_count = 0
    locus_len = 0
    name_regex = "(>.+)\n"
    name_regex_compile = re.compile(name_regex)

    for file_select in os.listdir(dir_of_aligns):
        #dict_of_loci_and_len_and_order[locus_count] = 
        full_path_to_file = dir_of_aligns + file_select
        open_file = open(full_path_to_file,'r')
        read_file = open_file.read()
        get_names = re.findall(name_regex_compile, read_file)
        if get_names:
            for name in get_names:
                if name not in dict_of_names_and_seqs:
                    dict_of_names_and_seqs[name] = []

                # TODO: HANDLE HAVING GAPS IN THE ALIGNMENT
                # CURRENTLY ONLY TAKES NUCLEOTIDES AND DEGENERATE BASES

                seq_grabber = name + "\n(\w|\n)+"
                seq_grabber_compile = re.compile(seq_grabber, re.S)
                seq_search = re.search(seq_grabber, read_file)
                split_name_and_seq = seq_search.group().split("\n", 1)
                contig_seq = split_name_and_seq[1].replace("\n", "")
                locus_len = len(contig_seq)
                dict_of_names_and_seqs[name].append(contig_seq)

        loci_info_dict = {}
        loci_info_dict[file_select] = locus_len
        dict_of_loci_and_len_and_order[locus_count] = loci_info_dict
        locus_count+=1

    output_file = open(args.out_file,'w')
    for name, seqs in dict_of_names_and_seqs.items():
        concat_seqs = ''.join(seqs)
        output_file.write(name)
        output_file.write("\n")
        output_file.write(concat_seqs)
        output_file.write("\n")
    output_file.close()
    
    print(dict_of_loci_and_len_and_order)
    #print(dict_of_names_and_seqs)



if __name__ == '__main__':
    main()
