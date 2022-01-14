#! /usr/bin/python3

import os
import argparse
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(prog='vcffixer test', \
        description='Automate downloading of high-throughput sequence data and updating of alignments using Extensiphy.')
    parser.add_argument('--ep_path')
    return parser.parse_args()

def main():
    args = parse_args()
    ep_path = args.ep_path

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]
    hand_checked_vcf = absolute_path + '/example_vcffixer.vcf' #can use full path here
    hand_checked_align = absolute_path + '/example_vcffixer.fa' #and here
    output_file = absolute_path + '/vcffixer_test_output.fa' #and here
    print(absolute_path)

    test_vcffixer(hand_checked_vcf, hand_checked_align, output_file, ep_path)

    check_test_output(hand_checked_align, output_file)

def test_vcffixer(vcf_, align_, output_, path):
    """Runs a test of vcffixer to ensure the output is as expected."""

    vfix = subprocess.Popen([path + "/modules/vcffixer.py", "--vcf_file", vcf_, "--align_file", align_, "--out_file", output_], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # print(vfix.communicate())

    assert os.path.exists(output_)

def check_test_output(in_align_, out_align_):

    input = open(in_align_, 'r')
    read_input = input.read()
    i_split = read_input.split("\n")

    output = open(out_align_, 'r')
    read_output = output.read()
    o_split = read_output.split("\n")

    print(len(i_split[1]))
    print(len(o_split[1]))
    assert (len(i_split[1]) + 2) == len(o_split[1])




if __name__ == '__main__':
    main()
