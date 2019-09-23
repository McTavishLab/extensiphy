#! /usr/bin/python

import os
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--msa_folder')
    parser.add_argument('--out_file')
    parser.add_argument('--suffix')
    return parser.parse_args()


def main():
    args = parse_args()

    msa_list = os.listdir(args.msa_folder)
    print(msa_list)

    name_finder = "(>.*\n)"
    name_finder_compiled = re.compile(name_finder)
   
    seq_dir = {}
    file_count = 0
    seq_count = 0
    for file in msa_list:
        if file.endswith(args.suffix):
            file_count+=1
            open_file = open(args.msa_folder + file, "r")
            read_file = open_file.read()
            name_findall = re.findall(name_finder, read_file)
            if file_count == 1:
                for name in name_findall:
                    newline_name_strip = name.strip("\n")
                    seq_dir[newline_name_strip] = ''
            #for name in name_findall:
                    seq_count+=1
                    seq_finder = name + "(.+)\n>"
                    seq_finder_compiled = re.compile(seq_finder, re.S)
                    seq_findall = re.findall(seq_finder_compiled, read_file)
                    name_check = name.strip("\n")
                    if name_check in seq_dir:
                        seq_dir[name_check] = '' + str(seq_findall)
                


            #if file_count == 1:
            #    for name in name_findall:
            #        newline_name_strip = name.strip("\n")
            #        seq_dir[newline_name_strip] = ''
                    

    print(seq_dir)
            









if __name__ == '__main__':
    main()
