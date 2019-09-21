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
    
    for file in msa_list:
        if file.endswith(args.suffix):
            open_file = open(args.msa_folder + file, "r")
            read_file = open_file.read()
            name_findall = re.findall(name_finder, read_file)
            print(name_findall)
            









if __name__ == '__main__':
    main()
