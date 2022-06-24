#! /usr/bin/env python3
import os
import re
import subprocess
import pathlib
from itertools import zip_longest
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

def mafft_align(cwd, unaligned_seqs_dir):
    """
    Use the mafft program to align the two sequences that have been identified as probably matching.
    """
    module_path = pathlib.Path(__file__).parent.absolute()

    output_dir = cwd + '/aligned_matched_seqs'
    if os.path.isdir(output_dir) == False:
        os.mkdir(output_dir)

    list_of_seqs = os.listdir(unaligned_seqs_dir)

    # mafft_command = ['mafft', '--op', '5', '--lexp', '-0.5', unaligned_seqs_dir + '/' +
    for num, seq_file in enumerate(list_of_seqs):
        input_file_path_and_name = unaligned_seqs_dir + '/' + seq_file
        output_file_path_and_name = output_dir + '/aligned_' + seq_file
        # mafft_command = ['mafft', '--op', '5', '--lexp', '-0.5', unaligned_seqs_dir + '/' + seq_file + ' > ' + output_dir + '/aligned_' + seq_file]
        mafft_command = ['mafft', '--op', '5', '--lexp', '-0.5', input_file_path_and_name, '>', output_file_path_and_name]
        print(mafft_command)
        subprocess.run([str(module_path) + '/mafft_command.sh', input_file_path_and_name, output_file_path_and_name])
        # align_result = subprocess.run(mafft_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # print(align_result.stdout)

    return output_dir


def alignment_fixer(read_file):
    split_seqs = read_file.split('>')

    output = []
    seq_1 = []
    seq_2 = []
    seq_count = 0
    for seq_and_name in split_seqs:
        if len(seq_and_name) > 0:
            seq_count+=1
            split_seq_and_name = seq_and_name.split('\n', 1)
            #print(split_seq_and_name)
            name = split_seq_and_name[0]
            seq = split_seq_and_name[1]

            fixed_seq = seq.replace('\n','')
            if seq_count == 1:
                seq_1.append(name)
                seq_1.append(fixed_seq)
            elif seq_count == 2:
                seq_2.append(name)
                seq_2.append(fixed_seq)

    output.append(seq_1)
    output.append(seq_2)

    return output
