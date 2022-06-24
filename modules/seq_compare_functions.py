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


def nuc_counter(sequence):
    nuc_count = 0
    for nuc in sequence[1]:
        if nuc != '-':
            nuc_count+=1

    return nuc_count

def trim_gaps(short_seq):
    output = []
    leading_gaps = 0
    trailing_gaps = 0

    leading_gap_match = "^(-*)\w"
    trailing_gap_match = "\w(-*)$"
    compile_leading_match = re.compile(leading_gap_match)
    compile_trailing_match = re.compile(trailing_gap_match)
    find_leading = re.findall(compile_leading_match, short_seq)
    find_trailing = re.findall(compile_trailing_match, short_seq)

    if compile_leading_match:
        #print("found leading")
        #print(find_leading)
        leading_gaps = len(find_leading[0])

        print("leading gap size", leading_gaps)

    if compile_trailing_match:
        #print("found trailing")
        #print(find_trailing)
        trailing_gaps = len(find_trailing[0])
        print("trailing gap size", trailing_gaps)

    output.append(leading_gaps)
    output.append(trailing_gaps)

    return output

def check_for_gaps(processed_seq, gaps, nucs):

    output = {}
    gap_positions = []
    nuc_positions = []
    for num, nuc in enumerate(processed_seq):
        if nuc.upper() in nucs:
            nuc_positions.append(num)
        elif nuc.upper() in gaps:
            gap_positions.append(num)
    contiguous_nucs = []

    num_gaps = len(gap_positions)
    num_nucs = len(nuc_positions)

    # check if the number of gaps in this sequence (hypothetically, the shorter sequence)
    # is greater or equal to 30% of the total number of nucleotides
    # if so, we've got a problem. Fix
    contiguous_nucs.append(nuc_positions[0])
    if num_gaps >= num_nucs * 0.3:
        for num, pos in enumerate(nuc_positions[1:]):
            num = num + 1
            if pos == nuc_positions[num-1] + 1:
                # print(pos)
                pass
            else:
                contiguous_nucs.append(nuc_positions[num-1])
                contiguous_nucs.append(nuc_positions[num])
        contiguous_nucs.append(max(nuc_positions))
        print(contiguous_nucs)
        largest_group = 0
        selected_group = []
        for num, item in enumerate(grouper(contiguous_nucs,2)):
            total_len = item[1] - item[0]
            if total_len > largest_group:
                largest_group = total_len
                selected_group.append([item[0],item[1]])
        print(selected_group[-1])
        output['anomalies'] = selected_group[-1]
        return output

    elif num_gaps < num_nucs * 0.3:
        output['no_anomalies'] = 0
        return output



def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)



def check_alignment(list_of_paired_nucs):
    total_nucs = 0
    identical_nucs = 0
    non_identical_nucs = 0
    gaps = 0
    gap_positions = []
    identical_positions = []
    non_identical_positions = []
    gap_set = ['-','N']
    nuc_set = ['A','C','G','T']
    output = []
    for num, pair in enumerate(list_of_paired_nucs):
        assert len(pair) == 2
        total_nucs+=1
        #print(pair[0])
        #print(pair[1])
        if str(pair[0].upper()) in gap_set or str(pair[1].upper()) in gap_set:
            #if pair[0] in gap_set:
            #    print(pair[0])
            #elif pair[1] in gap_set:
            #    print(pair[1])
            gaps+=1
            gap_positions.append(num)
            #print('gap')

        elif pair[0] and pair[1] not in gap_set:
            #print('no gap')
            if pair[0].upper() == pair[1].upper():
                identical_nucs+=1
                identical_positions.append(num)

            elif pair[0].upper() != pair[1].upper():
                non_identical_nucs+=1
                non_identical_positions.append(num)


    output.append(identical_nucs)
    output.append(non_identical_nucs)
    output.append(gaps)
    output.append(total_nucs)
    output.append(identical_positions)
    output.append(non_identical_positions)
    output.append(gap_positions)

    #return identical_nucs
    return output

def comparison(list_of_list_of_seqs):

    seq_1 = None
    seq_2 = None
    gap_set = ['-','N']
    nuc_set = ['A','C','G','T']
    seq_count = 0
    for list_of_seqs in list_of_list_of_seqs:
        seq_count+=1
        if seq_count == 1:
            seq_1 = list_of_seqs
        elif seq_count == 2:
            seq_2 = list_of_seqs

    #print(seq_1)
    #print(seq_2)

    count_nucs_1 = nuc_counter(seq_1)
    count_nucs_2 = nuc_counter(seq_2)

    # print(count_nucs_1)
    # print(count_nucs_2)

    nuc_lens = [count_nucs_1, count_nucs_2]
    small_seq = min(nuc_lens)

    shorter = None
    longer = None
    if small_seq == count_nucs_1:
        shorter = seq_1
        longer = seq_2

    elif small_seq == count_nucs_2:
        shorter = seq_2
        longer = seq_1

    # print("length of longer seq", len(shorter[1]))
    # print("length of shorter seq", len(longer[1]))

    shorter_seq = shorter[1]
    longer_seq = longer[1]

    len_short = nuc_counter(shorter)
    #len_long = nuc(longer)
    len_long = len(longer[1])
    #print(longer)

    get_gaps = trim_gaps(shorter_seq)
    #print(get_gaps)

    trimed_shorter = ''
    trimmed_longer = ''
    #print("waffle1")
    if get_gaps[1] != 0:

        #print(shorter_seq[get_gaps:-get_gaps[1]])
        trimmed_shorter = shorter_seq[get_gaps[0]:-get_gaps[1]]
        trimmed_longer = longer_seq[get_gaps[0]:-get_gaps[1]]
        #print(trimmed_shorter)
    elif get_gaps[1] == 0:
        #print(shorter_seq[get_gaps[0]:])
        trimmed_shorter = shorter_seq[get_gaps[0]:]
        trimmed_longer = longer_seq[get_gaps[0]:]
    #print("waffle2")

    gap_test = check_for_gaps(trimmed_shorter,gap_set, nuc_set)
    print(gap_test)
    for key, value in gap_test.items():
        if key == 'no_anomalies':
            print('no anomalous spread of nucleotides found.')
        elif key == 'anomalies':
            print('problem found')
            #TODO: add stuff here for handling of problem
            print(value[0])
            print(value[1])
            trimmed_shorter = trimmed_shorter[value[0]:value[1]]
            trimmed_longer = trimmed_longer[value[0]:value[1]]

    split_trimmed_short = list(trimmed_shorter)
    split_trimmed_long = list(trimmed_longer)

    combined_positions = list(map(list, zip(split_trimmed_long, split_trimmed_short)))
    analyze_alignment = check_alignment(combined_positions)

    #print(analyze_alignment)
    analyze_alignment.append(len_short)
    analyze_alignment.append(len_long)

    return analyze_alignment


def handle_comparison_outputs(cwd, comparison_info, output_file):
    print("identical nucs", comparison_info[0])
    print("non identical nucs", comparison_info[1])
    print("gaps", comparison_info[2])

#     output_file = open(output_file, 'w')
# #
#     # output_file.write(args.align_1)
#     output_file.write("\n")
#     output_file.write(">identical_nucs")
#     output_file.write("\n")
#     output_file.write(str(comparison_info[0]))
#     output_file.write("\n")
#     output_file.write(">non_identical_nucs")
#     output_file.write("\n")
#     output_file.write(str(comparison_info[1]))
#     output_file.write("\n")
#     output_file.write(">gaps")
#     output_file.write("\n")
#     output_file.write(str(comparison_info[2]))
#     output_file.write("\n")
#     output_file.write(">total_nucleotides")
#     output_file.write("\n")
#     output_file.write(str(comparison_info[3]))
#     output_file.write("\n")
#     output_file.write(">identical_nucs_positions")
#     output_file.write("\n")
#     output_file.write(str(comparison_info[4]))
#     output_file.write("\n")
#     output_file.write(">non_identical_nucs_positions")
#     output_file.write("\n")
#     output_file.write(str(comparison_info[5]))
#     output_file.write("\n")
#     output_file.write(">gaps_positions")
#     output_file.write("\n")
#     output_file.write(str(comparison_info[6]))
#     output_file.write("\n")
#     output_file.write(">unadjusted_short_seq_length")
#     output_file.write("\n")
#     output_file.write(str(comparison_info[7]))
#     output_file.write("\n")
#     output_file.write(">unadjusted_long_seq_length")
#     output_file.write("\n")
#     output_file.write(str(comparison_info[8]))
