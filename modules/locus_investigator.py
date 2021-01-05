#! /usr/bin/env python3
import argparse
import os
import re
import csv
# from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--locus')
    parser.add_argument('--position_file')
    return parser.parse_args()

def main():
    args = parse_args()

    input_file = os.path.realpath(args.position_file)

    with open(input_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        nucs_to_specific_cluster = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                nucs_to_specific_cluster = nucs_to_specific_cluster + (int(row[2]))
                if args.locus in row[1]:
                    print(row)
                    print(nucs_to_specific_cluster - int(row[2]))
                # print(f'\t{row[0]} {row[1]} {row[2]}.')
                line_count += 1
        print(f'Processed {line_count} lines.')

    



    
if __name__ == '__main__':
    main()