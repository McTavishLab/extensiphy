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

    name_regex = "(>.+)\n"
    name_regex_compile = re.compile(name_regex)

    for file_select in os.listdir(dir_of_aligns):
        full_path_to_file = dir_of_aligns + file_select
        open_file = open(full_path_to_file,'r')
        read_file = open_file.read()
        get_names = re.findall(name_regex_compile, read_file)
        if get_names:
            for name in get_names:
                if name not in dict_of_names_and_seqs:
                    dict_of_names_and_seqs[name] = []

    print(dict_of_names_and_seqs)

    # TODO NOW START ADDING SEQUENCES TO THE LISTS ATTACHED TO THE NAMES IN THE DICT
    # ADDITIONALLY, KEEP TRACK OF LENGTHS AND OUTPUT LENGTHS AND ORDER OF SEQUENCES AS THEY ARE ADDED TO THE FINAL MSA
        



if __name__ == '__main__':
    main()
