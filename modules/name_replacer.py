#! /usr/bin/env python3
import argparse
import os
from os import walk
import re
import pandas as pd
import csv

def process_reads(input_file, data_container, name_file_input, tail_1, tail_2):
    replacement_name_count = 0

    for (dirpath, dirnames, filenames) in walk(input_file):
        print(filenames)
        for file in filenames:
            file_1 = ''
            file_2 = ''
            if tail_1 in file:
                file_1 = file

                replacement_name_count+=1
                
                for potential_second_file in filenames:
                    if potential_second_file == file_1.replace(tail_1, tail_2):
                        file_2 = potential_second_file
                        print(dirpath + '/' + file_1)
                        print(dirpath + '/' + file_2)
                        print("next pair")
                        data_container.append([replacement_name_count, file_1, file_2, tail_1, tail_2])
                        
                        replacement_file_1 = str(replacement_name_count) + '_' + tail_1
                        replacement_file_2 = str(replacement_name_count) + '_' + tail_2

                        os.rename(dirpath + '/' + file_1, replacement_file_1)
                        os.rename(dirpath + '/' + file_2, replacement_file_2)

    # return data_container

def write_csv(df, output_file_name):
    if output_file_name == "NONE":
        df.to_csv('original_read_file_names.csv', index=False)


def parse_args():
    parser = argparse.ArgumentParser()
    # parser.add_argument('-r', action='store_true', help='parses fasta files and creates a name replacement for each file name that is then used in rapup')
    # parser.add_argument('-a', action='store_true', help='parses fasta files and creates a name replacement for each file name that is then used in rapup')
    # parser.add_argument('-ra', action='store_true', help='parses fasta files and creates a name replacement for each file name that is then used in rapup')
    parser.add_argument('--align', type=str, default="NONE")
    parser.add_argument('--read_dir', type=str, default="NONE")
    parser.add_argument('--tail_1', type=str, default="NONE")
    parser.add_argument('--tail_2', type=str, default="NONE")
    parser.add_argument('--name_file', type=str, default="NONE")
    parser.add_argument('--outfile', type=str, default="NONE")
    return parser.parse_args()

def main():
    args = parse_args()

    columns = ['new_name', 'old_file_1', 'old_file_2', 'tail_1', 'tail_2']
    # df_ = pd.DataFrame(index=index, columns=columns)

    data = []
    
    if args.read_dir != "NONE" and args.name_file == "NONE":
        process_reads(args.read_dir, data, "NONE", args.tail_1, args.tail_2)

    df_ = pd.DataFrame(data, columns=columns)
    print(df_)

    write_csv(df_, "NONE")

if __name__ == '__main__':
    main()