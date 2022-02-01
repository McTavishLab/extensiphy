#! /usr/bin/python3

import os
import argparse
import subprocess
import time
import pytest

def parse_args():
    parser = argparse.ArgumentParser(prog='vcffixer test', \
        description='This program tests the vcffixer.py module of Extensiphy. \
        Testing currently ensures that the length of the sequence output by vcffixer.py \
        matches the length stated in the vcf used by Extensiphy. \
        EXAMPLE COMMAND: vcffixer_test.py \
        EXAMPLE PYTEST COMMAND: pytest vcffixer_test.py')
    # parser.add_argument('--ep_path', help='Absolute path to your Extensiphy directory.')
    return parser.parse_args()

def main():
    args = parse_args()
    # ep_path = args.ep_path

    # split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    # absolute_path = split_path_and_name[0]
    # ep_path = '/' + absolute_path.strip('/tests')
    # hand_checked_vcf = absolute_path + '/example_vcffixer.vcf' #can use full path here
    # hand_checked_align = absolute_path + '/example_vcffixer.fa' #and here
    # output_file = absolute_path + '/TEST_OUTPUT_vcffixer.fa' #and here
    # print(absolute_path)
    # print(ep_path)

    # test_vcffixer(hand_checked_vcf, hand_checked_align, output_file, ep_path)
    test_vcffixer()

    # check_test_output(hand_checked_align, output_file)
    test_output()

#vcf_, align_, output_, path
def test_vcffixer():
    """Runs a test of vcffixer to ensure the output is as expected."""
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]
    path = '/' + absolute_path.strip('/tests')
    vcf = absolute_path + '/example_vcffixer.vcf' #can use full path here
    align = absolute_path + '/example_vcffixer.fa' #and here
    output = absolute_path + '/TEST_OUTPUT_vcffixer.fa'

    print(path + "/modules/vcffixer.py", "--vcf_file", vcf, "--align_file", align, "--out_file", output)
    vfix = subprocess.Popen([path + "/modules/vcffixer.py", "--vcf_file", vcf, "--align_file", align, "--out_file", output], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # print(vfix.communicate())

    print("vcffixer.py run complete. Checkout output.")

    # This sleep line helps to avoice Race Conditions of the output file
    time.sleep(5)
    assert os.path.exists(output)

def test_output():

    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]
    align = absolute_path + '/example_vcffixer.fa' #and here
    output = absolute_path + '/TEST_OUTPUT_vcffixer.fa'

    input = open(align, 'r')
    read_input = input.read()
    i_split = read_input.split("\n")

    output = open(output, 'r')
    read_output = output.read()
    o_split = read_output.split("\n")

    print(len(i_split[1]))
    print(len(o_split[1]))
    assert (len(i_split[1]) + 2) == len(o_split[1])

    print("vcffixer.py test: PASSED.")


if __name__ == '__main__':
    main()
