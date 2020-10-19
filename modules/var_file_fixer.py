#! /usr/bin/python

import os
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--var_file')
    parser.add_argument('--out_file')
    return parser.parse_args()

def fix_var_file(file_contents):
    
    regex = '_=posix(.+)'
    compile_regex = re.compile(regex)

    replaced_newlines = file_contents.replace('\n', '$')
    # print(replaced_newlines)

    find_non_standard_vars = re.findall(compile_regex, replaced_newlines)
    if find_non_standard_vars:
        # print(find_non_standard_vars)
        output_vars = find_non_standard_vars[0].replace('$','\n')
        # print(output_vars)
        assert len(output_vars) > 1
        return output_vars


def main():
    args = parse_args()

    input = open(args.var_file,'r')
    read_input = input.read()
    # print(read_input)

    fix = fix_var_file(read_input)

    output_file = open(args.out_file, 'w')
    output_file.write(fix)
    output_file.close()

    input.close()

if __name__ == '__main__':
    main()