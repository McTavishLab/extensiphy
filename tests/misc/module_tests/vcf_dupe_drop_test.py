#! /usr/bin/python3

import os
import argparse
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(prog='vcffixer test', \
        description='This program tests the vcf_duplicate_dropper.py module of Extensiphy. \
        Testing currently ensures that the length of the sequence output by vcffixer.py \
        matches the length stated in the vcf used by Extensiphy. \
        EXAMPLE COMMAND: vcf_duplicate_dropper.py --ep_path [path to extensiphy]')
    parser.add_argument('--ep_path', help='Absolute path to your Extensiphy directory.')
    return parser.parse_args()

def main():
    args = parse_args()

    ep_path = args.ep_path
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]


if __name__ == '__main__':
    main()
