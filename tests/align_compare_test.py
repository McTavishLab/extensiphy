#! /usr/bin/python3

import os
import argparse
import subprocess
import time
import pytest
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(prog='align_compare test', \
        description='This program tests the align_compare.py module of Extensiphy. \
        Test currently checks that the output dataframe is exactly as expected \
        EXAMPLE COMMAND: align_compare_test.py \
        EXAMPLE PYTEST COMMAND: pytest align_compare_test.py')
    # parser.add_argument('--ep_path', help='Absolute path to your Extensiphy directory.')
    return parser.parse_args()

def main():
    args = parse_args()

    test_align_compare()

    test_output()



def test_align_compare():
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]
    path = '/' + absolute_path.strip('/tests')
    inputs_path = path + '/testdata'
    align_1 = inputs_path + '/combo.fas' #can use full path here
    align_2 = inputs_path + '/alt.fas' #and here
    output = absolute_path + '/TEST_OUTPUT_align_compare.csv'

    compare = subprocess.Popen([path + "/modules/align_compare.py", "--align_1", align_1, "--align_2", align_2, "--out_file", output], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    print("vcffixer.py run complete. Checkout output.")

    # This sleep line helps to avoice Race Conditions of the output file
    time.sleep(5)
    assert os.path.exists(output)

def test_output():
    split_path_and_name = os.path.realpath(__file__).rsplit('/',1)
    absolute_path = split_path_and_name[0]
    output = absolute_path + '/TEST_OUTPUT_align_compare.csv'

    hand_checked_list = [['taxon', 'position', 'file_one_base', 'file_two_base'], ['taxon_25', 37, 'T', 'C'], ['taxon_25', 172, 'G', 'T'], ['taxon_24', 10137, 'T', 'G']]

    hand_checked_df = pd.DataFrame(hand_checked_list[1:], columns=hand_checked_list[0])

    test_df = pd.read_csv(output, delimiter=',', index_col=0)

    print(hand_checked_df)
    print(test_df)

    assert hand_checked_df.equals(test_df)


if __name__ == '__main__':
    main()
