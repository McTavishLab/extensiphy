#!/usr/bin/env python3
#Functions for manipulating fasta files in various ways

import os
import argparse
import re
import random
from Bio.Seq import Seq

def split_fasta_into_seqs(align_file, out_dir, out_file_locus_name):
    """
    Function that splits a multiple sequence alignment file in fasta format \
    into individual sequences in individual files.
    """
    with open(align_file, 'r') as raw_seqs:

        read_seq = raw_seqs.read()
        split_seqs = read_seq.split(">")
        #print(split_seqs[1])
        for name_and_seq in split_seqs:
            if len(name_and_seq) > 2:
                split_name = name_and_seq.split("\n", 1)
                #print(split_name)
                name = split_name[0]
                #print(name)
                seq = split_name[1]
                cluster_name = out_file_locus_name
                full_fasta_name = out_dir + '/' + cluster_name + "_" + name
                open_file = open(full_fasta_name, "w")
                open_file.write(">")
                open_file.write(name_and_seq)


def split_multiple_fastas_into_seqs(cwd, align_dir, align_suffix):
    """
    Uses the split_fasta_into_seqs function to split up many fasta files into individual sequences. \
    Ideal for use when analyzing the bases of a concatenated alignment updated with Extensiphy.
    """

    if os.path.isdir(cwd + '/split_locus_files') == False:
        os.mkdir(cwd + '/split_locus_files')

    # os.chdir(cwd + '/split_locus_files')

    out_dir = cwd + '/split_locus_files'

    path = os.path.abspath(align_dir)

    list_of_files = os.listdir(path)

    for file in list_of_files:
        path_to_align = path + '/' + file
        if os.path.isfile(path_to_align):
            if file.endswith(align_suffix):
                suffix_removed_align_name = file.strip(align_suffix)
                # print(suffix_removed_align_name)

                split_fasta_into_seqs(path_to_align, out_dir, suffix_removed_align_name)
            else:
                print("FILE NOT ENDING IN SUFFIX: ", file)

        else:
            print("DIRECTORY: ", path_to_align)


def make_genomes_contiguous(cwd, genome_dir, suffix):
    """
    de novo assembled genomes are often split into "chromosomes" or other contigs.
    To simplify the search for specific loci, this function concatenates the genome into a single sequences.
    """
    list_of_genomes = os.listdir(genome_dir)
    abs_path_to_genomes = os.path.abspath(genome_dir)

    if os.path.isdir('contiguous_genome_files') == False:
        os.mkdir(cwd + '/contiguous_genome_files')

    os.chdir(cwd + '/contiguous_genome_files')

    for genome_file in list_of_genomes:
        if os.path.isfile:
            if genome_file.endswith(suffix):

                contiguous_seq_file = []
                new_seq = []

                # print(genome_file)
                stripped_suffix_genome = genome_file.strip(suffix)
                print(stripped_suffix_genome)

                contiguous_seq_file.append('>' + stripped_suffix_genome + '\n')

                file_contents = open(abs_path_to_genomes + '/' + genome_file, 'r').readlines()
                for line in file_contents:
                    if not line.startswith('>'):
                        new_seq.append(line.strip('\n'))

                #Join the sequences together to make a contiguos sequence
                contiguous_seq = ''.join(new_seq)
                # print(contiguous_seq)

                contiguous_seq_file.append(contiguous_seq)

                #Write new seq to file
                output = open('contiguous_' + genome_file, 'w')
                for item in contiguous_seq_file:
                    output.write(item)
                output.write('\n')
                output.close()

    os.chdir(cwd)


def combine_genomes_with_loci(loci_dir, genome_dir):
    """
    Identify the shared names of the genomes and the sequences in the cluster directory. \
    Combine the appropriate sequences so that the loci can be aligned to the genomes.
    """

    list_of_genomes = os.listdir(genome_dir)
    abs_path_to_genomes = os.path.abspath(genome_dir)

    list_of_loci = os.listdir(loci_dir)
    abs_path_to_loci = os.path.abspath(loci_dir)

    for genome_file in list_of_genomes:
        print(genome_file)


def seq_converter(seq):
    print("seq_converter")

    output = {}

    seq_split = seq.split('\n', 1)
    label = seq_split[0]
    seq = seq_split[1]

    len_seq = len(seq_split[1])
    shorter = seq.replace('\n','')

    main_short = Seq(shorter)
    short_comp = main_short.complement()
    short_reverse = main_short[::-1]
    short_rev_comp = main_short.reverse_complement()

    assert len(main_short) > 0
    assert len(short_comp) > 0
    assert len(short_reverse) > 0
    assert len(short_rev_comp) > 0

    assert len(main_short) == len(short_comp) == len(short_reverse) == len(short_rev_comp)

    output['original'] = str(main_short)
    output['complement'] = str(short_comp)
    output['reverse'] = str(short_reverse)
    output['reverse_complement'] = str(short_rev_comp)

    return output

def generate_random_kmer_positions(seq_len):
    print("generate_random_kmer_positions")
    output = []

    for i in range(0,50):
        kmer_start = random.randint(0,seq_len)
        if(kmer_start + 14 <= seq_len):
            output.append(kmer_start)
    return output

# def generate_kmer_positions(seq_len):
#     output = []

#     for i in range(0,50):
#         kmer_start = random.randint(0,seq_len)
#         if(kmer_start + 14 <= seq_len):
#             output.append(kmer_start)
#     return output

def find_match(long_seq, dict_of_seqs):
    orig_matches = 0
    comp_matches = 0
    rev_matches = 0
    rev_comp_matches = 0
    compare_nums = []


    for orientation, seq in dict_of_seqs.items():
        seq_len = len(seq)

        get_kmer_starts = generate_random_kmer_positions(seq_len)

        for position in get_kmer_starts:
            kmer = seq[position:position + 14]
            compile_kmer = re.compile(kmer.upper())
            find_kmer = re.findall(compile_kmer, long_seq.upper())
            if find_kmer:
                if orientation == 'original':
                    orig_matches+=1
                elif orientation == 'complement':
                    comp_matches+=1
                elif orientation == 'reverse':
                    rev_matches+=1
                elif orientation == 'reverse_complement':
                    rev_comp_matches+=1

    compare_nums.append(orig_matches)
    compare_nums.append(comp_matches)
    compare_nums.append(rev_matches)
    compare_nums.append(rev_comp_matches)

    print(compare_nums)

    if max(compare_nums) == orig_matches:
        return 'original'
    elif max(compare_nums) == comp_matches:
        return 'complement'
    elif max(compare_nums) == rev_matches:
        return 'reverse'
    elif max(compare_nums) == rev_comp_matches:
        return 'reverse_complement'


def find_boundaries(manipulated_seq, long_seq):
    print("find_boundaries")
    output_kmers = []
    front_kmers = []
    back_kmers = []
    kmer_size = 50
    num_kmers = 20
    current_kmer_position = 0
    end_kmer_position = -1

    print(len(manipulated_seq))
    print('BEGINNING SECTION')
    for num in range(0,num_kmers):
        kmer = manipulated_seq[current_kmer_position:current_kmer_position + kmer_size]
        # print(kmer)
        # front_kmers.append(kmer)
        current_kmer_position = current_kmer_position + kmer_size
        compile_kmer = re.compile(kmer.upper())
        find_kmer = re.search(compile_kmer, long_seq.upper())
        if find_kmer:
            find_kmer = find_kmer.start()
            # print(find_kmer)
            front_kmers.append(find_kmer)
        else:
            front_kmers.append(0)

    print("END SECTION ")

    for num in range(0,num_kmers):
        kmer = manipulated_seq[end_kmer_position - kmer_size : end_kmer_position]
        # print(kmer)
        # back_kmers.append(kmer)
        end_kmer_position = end_kmer_position - kmer_size
        compile_kmer = re.compile(kmer.upper())
        find_kmer = re.search(compile_kmer, long_seq.upper())
        if find_kmer:
            find_kmer = find_kmer.start()
            # print(find_kmer)
            back_kmers.append(find_kmer)
        else:
            # back_kmers.append('-')
            back_kmers.append(0)



    print(front_kmers)
    # print(manipulated_seq)
    print(back_kmers)
    assert len(front_kmers) == num_kmers
    assert len(back_kmers) == num_kmers

    refine_start = refine_potential_boundaries(front_kmers, kmer_size, "start")
    refine_stop = refine_potential_boundaries(back_kmers, kmer_size, "stop")

    print(refine_start, refine_stop)
    print("WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW")

    output_kmers.append(refine_start)
    output_kmers.append(refine_stop)
    return output_kmers

def refine_potential_boundaries(kmer_match_list, kmer_len, orientation):
    print("refine_potential_boundaries")
    print("REFINE BOUNDARIES")
    regions = {}
    region_count = 0
    kmer_count = 0
    current_region = []
    longest_region = 0
    longest_region_value = 0
    start = 0
    output_position = 0
    # while kmer_count != len(kmer_match_list) - 1:
    for num, kmer in enumerate(kmer_match_list):
        if num + 1 < len(kmer_match_list):
            # print(kmer)
            # print(kmer_match_list[num + 1])
            # print("WAAAAAAAAAAAAAAAAAAAAAAAAAFFLE")

            if kmer_match_list[kmer_count] + kmer_len == kmer_match_list[kmer_count + 1] or kmer_match_list[kmer_count] - kmer_len == kmer_match_list[kmer_count + 1]:
                current_region.append(kmer_match_list[kmer_count])
                # regions[region_count] = []
            elif kmer_match_list[kmer_count] + kmer_len != kmer_match_list[kmer_count + 1] or kmer_match_list[kmer_count] - kmer_len == kmer_match_list[kmer_count + 1]:
                current_region.append(kmer_match_list[kmer_count])
                regions[region_count] = current_region
                current_region = []
                region_count+=1

            kmer_count+=1
    regions[region_count] = current_region
    # for key, value in regions.items():

    #     if value[-1] + kmer_len == kmer_match_list[-1]:
    #         value.append(kmer_match_list[-1])

    print(regions)
    for key, value in regions.items():
        if len(value) > longest_region_value:
            longest_region_value = len(value)
            longest_region = key

    longest_match_region = regions.get(longest_region)

    if orientation == "start":
        start = min(longest_match_region)

    elif orientation == "stop":
        start = max(longest_match_region)

    for num, position in enumerate(kmer_match_list):
        if position == start:
            if orientation == "start":
                output_position = start - ((num + 1) * kmer_len)
            elif orientation == "stop":
                output_position = start + ((num + 1 )* kmer_len)
    print(output_position)
    return output_position
    # print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")


def trim_boundaries(kmer_lists, long_seq, kmer_len):
    print("trim_boundaries")
    contiguous_front_positions = []
    contiguous_back_positions = []
    seq_front_kmers = kmer_lists[0]
    seq_back_kmers = kmer_lists[1]
    buffered_start_position = 0
    buffered_stop_position = 0
    step_count = 1
    buffer_size = 150
    for num, kmer_start in enumerate(seq_front_kmers):
        # print(num)
        # print(kmer_start)
        if num == 0:
            continue
        elif num > 0:
            if type(kmer_start) == int and type(seq_front_kmers[num - 1]) == int and kmer_start == seq_front_kmers[num - 1] + kmer_len:
                contiguous_front_positions.append(num)
    # print(contiguous_front_positions)
    if len(contiguous_front_positions) == 0:
        contiguous_front_positions.append('no_useful_matches')
    else:
        earliest_starting_position = min(contiguous_front_positions)

        # calculate the number of steps between the contiguous starting point and the actual start of the kmers
        step_number = earliest_starting_position - step_count

        starting_seqence_position = seq_front_kmers[min(contiguous_front_positions)] - (kmer_len * step_number)
        buffered_start_position = starting_seqence_position - buffer_size
        print("starting position: ", buffered_start_position)

    # CALCULATE ENDING REGION AND BUFFER
    print("CALCULATE ENDING REGION AND BUFFER")
    for num, kmer_start in enumerate(seq_back_kmers):
        # print(num)
        # print(kmer_start)
        if num == 0:
            continue
        elif num > 0:
            if type(kmer_start) == int and type(seq_front_kmers[num - 1]) == int and kmer_start == seq_back_kmers[num - 1] - kmer_len:
                contiguous_back_positions.append(num)

    # print(contiguous_back_positions)
    if len(contiguous_back_positions) == 0:
        contiguous_back_positions.append('no_useful_matches')
    else:
        # print("calculating backside")
        earliest_stopping_position = min(contiguous_back_positions)
        # print(earliest_stopping_position)

        # calculate the number of steps between the contiguous starting point and the actual start of the kmers
        step_number = earliest_stopping_position - step_count
        # print(step_number)

        stopping_seqence_position = seq_back_kmers[min(contiguous_back_positions)] + (kmer_len * step_number)
        # print(stopping_seqence_position)

        buffered_stop_position = stopping_seqence_position + buffer_size
        print("stopping position: ", buffered_stop_position)

    if contiguous_front_positions[0] != 'no_useful_matches' and contiguous_back_positions[0] != 'no_useful_matches':
        print("LENGTH OF SEQUENCE REGION: ", buffered_stop_position - buffered_start_position)
        output = [buffered_start_position, buffered_stop_position]
        return output
    elif contiguous_front_positions[0] == 'no_useful_matches':
        output = ['no_useful_matches', buffered_stop_position]
        return output
    elif contiguous_back_positions[0] == 'no_useful_matches':
        output = [buffered_start_position, 'no_useful_matches']
        return output


def long_seq_trimmer(long_seq, len_short_seq, start_stop_positions_list):
    print("long_seq_trimmer")
    # start = start_stop_positions_list[0]
    # stop = start_stop_positions_list[1]
    start = None
    stop = None

    # if type(start) == int and type(stop) == int:
    if type(start_stop_positions_list[0]) == int and type(start_stop_positions_list[1]) == int:
        start = min(start_stop_positions_list)
        stop = max(start_stop_positions_list)
        output = long_seq[start : stop]
        # print(output)
        return output
    # elif start == 'no_useful_matches':
    if start_stop_positions_list[0] == 'no_useful_matches':
        stop = start_stop_positions_list[1]
        start = stop - (len_short_seq + 500)
        output = long_seq[start : stop]
        # print(output)
        return output
    # elif stop == 'no_useful_matches':
    if start_stop_positions_list[1] == 'no_useful_matches':
        start = start_stop_positions_list[0]
        stop = start + (len_short_seq + 500)
        output = long_seq[start : stop]
        # print(output)
        return output

# Test to assert length of region calculated as locus boundaries
# is similar in size to the length of the short locus
def assess_boundaries(short_seq, boundaries):
    print("assess_boundaries")
    split_seq_and_name = short_seq.split('\n', 1)
    seq = split_seq_and_name[1]
    short_len = len(seq)
    boundary_len = boundaries[1] - boundaries[0]
    assert boundary_len >= .999 * short_len


def match_long_with_loci(manip_seq_path, long_seq_path, output_dir):
    print("match_long_with_loci")
    kmer_len = 50
    manip_folder_contents = os.listdir(manip_seq_path)
    long_seqs_folder_contents = os.listdir(long_seq_path)
    file_info_regex = r'(cluster\d+)-_(.+)$'
    file_info_compile = re.compile(file_info_regex)
    long_name_regex = r'(\w+)\.fasta$'
    long_name_compile = re.compile(long_name_regex)

    num_long_files = len(long_seqs_folder_contents)
    num_short_files = len(manip_folder_contents)
    # print("number of files")
    # print(num_long_files)
    # print(num_short_files)

    manip_file_count = 0
    long_file_count = 0

    print("iterate over manip files")
    for manip_file in manip_folder_contents:
        # print(manip_file)
        find_info = re.findall(file_info_compile, manip_file)
        if find_info:
            print(find_info)
    #         manip_file_count+=1
    #         manip_taxon = find_info[0][1]
    #         manip_locus = find_info[0][0]
    #         # print(manip_taxon)
    #         # print(manip_locus)
    #         # print("iterate over long files")
    #         long_file_count = 0
    #         for long_seq in long_seqs_folder_contents:
    #             find_long_info = re.findall(long_name_compile, long_seq)
    #             #print("finding long seq info")
    #             if find_long_info:
    #                 long_file_count+=1
    #                 long_seq_name = find_long_info[0]
    #                 if long_seq_name == manip_taxon:
    #                     print("taxon match")
    #                     # print(long_seq)
    #                     # print(manip_taxon)
    #                     open_manip_file = open(manip_seq_path +'/'+ manip_file,'r')
    #                     read_manip_file = open_manip_file.read()
    #
    #                     open_long_seq = open(long_seq_path +'/'+ long_seq, 'r')
    #                     read_long_seq = open_long_seq.read()
    #
    #                     # print("files opened")
    #                     convert_manip = seq_converter(read_manip_file)
    #
    #                     long_seq_split = read_long_seq.split('\n', 1)
    #                     label = long_seq_split[0]
    #                     seq = long_seq_split[1]
    #
    #                     len_seq = len(long_seq_split[1])
    #                     long_contiguous = seq.replace('\n','')
    #
    #                     match_maker = find_match(long_contiguous, convert_manip)
    #                     print(match_maker)
    #
    #                     find_seq_location = find_boundaries(convert_manip[match_maker], long_contiguous)
    #
    #                     #check_boundaries = assess_boundaries(read_manip_file, find_seq_location)
    #
    #                     trimmed_long = long_contiguous[find_seq_location[0]:find_seq_location[1] + 100]
    #
    #                     # trim_long = trim_boundaries(find_seq_location, long_contiguous, kmer_len)
    #                     # print(trim_long)
    #
    #                     # produce_trimmed_seq = long_seq_trimmer(long_contiguous, len(convert_manip[match_maker]), trim_long)
    #
    #                     output_file_name = output_dir +'/'+ 'combined-' + manip_locus + '--' + manip_taxon
    #                     print(output_file_name)
    #                     output = open(output_file_name, 'w')
    #                     output.write(label)
    #                     output.write('\n')
    #                     output.write(trimmed_long)
    #                     output.write('\n')
    #                     output.write(label)
    #                     output.write('\n')
    #                     output.write(convert_manip[match_maker])
    #
    #                     output.close()
    #                     open_long_seq.close()
    #                     open_manip_file.close()
    #
    #                 # else:
    #                 #     #print("didnt find match between method sequences")
    #                 #     #print(manip_taxon + '    ' + long_seq_name)
    #                 #     long_file_count+=1
    #                 #     if long_file_count == num_long_files:
    #                 #         print("couldnt find match for this short file")
    #                 #         print(manip_taxon)
    #
    #
    #             else:
    #                 print("didnt find long match")
    #                 print(find_long_info)
    #
    #
    #     else:
    #         print("didnt find manip match")
    #         print(file_info_compile)
    #
    # print("number of files identified by regex methods")
    # print(manip_file_count)
    # print(long_file_count)
