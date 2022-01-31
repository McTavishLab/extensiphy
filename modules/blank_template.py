#!/usr/bin/env python3
# Blank python template using standard functional programming techniniques.
import os
import subprocess
import argparse

import sys
if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")

_DEBUG = 0
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_file')
    parser.add_argument('--align_file') ##consensus sequence
    parser.add_argument('--out_file')
    return parser.parse_args()


def main():
    args = parse_args()




if __name__ == '__main__':
    main()
