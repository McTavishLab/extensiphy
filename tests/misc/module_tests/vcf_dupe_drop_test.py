#! /usr/bin/python3

import os
import argparse
import subprocess
import time

def parse_args():
    parser = argparse.ArgumentParser(prog='vcf_duplicate_dropper test', \
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
    hand_checked_vcf = absolute_path + '/example_dupe_drop.vcf' #can use full path here
    output_file = absolute_path + '/vcf_drop_test_output.vcf' #and here

    test_dupe_dropper(ep_path, hand_checked_vcf, output_file)

    check_output(output_file)

    print("vcf_duplicate_dropper.py functioning normally.")

def test_dupe_dropper(path_, vcf_file_, output_):
    """Runs a test of vcf_duplicate_dropper to ensure the output is as expected."""

    vfix = subprocess.Popen([path_ + "/modules/vcf_duplicate_dropper.py", "--vcf_file", vcf_file_, "--out_file", output_], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # print(vfix.communicate())
    print("Test vcf produced.")

    #line helps prevent race condition
    time.sleep(5)
    assert os.path.exists(output_)

def check_output(output_):
    """Check the output and verify that no duplicate positions are found in the output"""
    output_pos_list = []
    # open_vcf = open(output_,'r')
    # read_vcf = open_vcf.readlines()
    line_count = 0
    # for line in read_vcf:
    #     if not line.startswith("#"):
    with open(output_) as read_vcf:
        for line in read_vcf:
            if not line.startswith("#"):
                line_count+=1
                # print(line)
                splitter = line.split()
                # print(splitter)

                position = splitter[1]
                if position in output_pos_list:
                    # print(output_pos_list)
                    print("ERROR. Duplicate position found after running vcf_duplicate_dropper.py")
                    print(line)
                    print("File line: ", line_count)
                    print("Exiting.")
                    exit()
                elif position not in output_pos_list:
                    output_pos_list.append(position)







if __name__ == '__main__':
    main()
