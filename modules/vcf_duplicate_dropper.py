#! /usr/bin/python

import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_file')
    parser.add_argument('--out_file')
    return parser.parse_args()


def main():
    args = parse_args()

    new_vcf = []
    duplicate_position_count = 0
    duplicate_position_check = 0
    with open(args.vcf_file) as vcf:
        for line_num, line in enumerate(vcf):

            if not line.startswith("#"):

                splitter = line.split()
                ref_name = splitter[0:1]
                pos = splitter[1:2]
                if duplicate_position_check != pos:
                    duplicate_position_check = pos
                    new_vcf.append(line)
                elif duplicate_position_check == pos:
                    print(line)
                # ref = splitter[3:4]
                # alt = splitter[4:5]
            elif line.startswith("#"):
                new_vcf.append(line)
    new_vcf_file = open(args.out_file, 'a+')
    for line in new_vcf:
        new_vcf_file.write(line)

    new_vcf_file.close()

if __name__ == '__main__':
    main()
