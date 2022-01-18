#! /usr/bin/python
# Program that dropes duplicate position entries from a VCF file.
#

import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_file')
    parser.add_argument('--out_file')
    return parser.parse_args()


def main():
    args = parse_args()

    drop_duplicate_positions(args.vcf_file, args.out_file)


def drop_duplicate_positions(vcf_file_, out_file_name_):
    """Function reads in VCF file, loops over the lines and notes positions. \
    If a position comes up more than once, only the first instance of that positon \
    is output to a new VCF file."""
    new_vcf = []
    duplicate_position_count = 0
    duplicate_position_check = 0
    with open(vcf_file_) as vcf:
        for line_num, line in enumerate(vcf):

            if not line.startswith("#"):

                splitter = line.split()
                # ref_name = splitter[0:1]
                ref_name = splitter[0]
                # pos = splitter[1:2]
                pos = int(splitter[1])
                if duplicate_position_check != pos:
                    duplicate_position_check = pos
                    new_vcf.append(line)
                elif duplicate_position_check == pos:
                    print(line)
                # ref = splitter[3:4]
                # alt = splitter[4:5]
            elif line.startswith("#"):
                new_vcf.append(line)
    # new_vcf_file = open(out_file_name_, 'a+')
    new_vcf_file = open(out_file_name_, 'w')
    for line in new_vcf:
        new_vcf_file.write(line)

    new_vcf_file.close()

if __name__ == '__main__':
    main()
