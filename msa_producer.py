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

    for file in os.listdir(dir_of_aligns):
        print(file)



if __name__ == '__main__':
    main()
