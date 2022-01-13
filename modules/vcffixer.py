#!/usr/bin/env python
# Program handles an issue found in VCF files that use the reference nucleotide
# in the flattened consensus sequence when the reference is an 'N'
# instead of the nucleotide found in the reads
# USE: vcffixer.py --vcf_file [VCF_FILE] --align_file [ALIGNMENT_FILE]
# TODO: Write this program in bash to reduce processing time.
# TODO: Contiguity of the VCF isnt necessary, just check the Ns to make sure there isnt an alt nucleotide at that position
import os
import subprocess
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_file')
    parser.add_argument('--align_file')
    parser.add_argument('--out_file')
    return parser.parse_args()


def main():
    args = parse_args()

    nuc_set = ['A', 'C', 'G', 'T']

    # create an empty list that will hold the fixed sequence
    ref_sequence_list = []
    ref_sequence_dict = {}
    output_list = []

    # initialize dictionary for location of nucleotide to be replaced
    # and the nucleotide that will go there
    replace_nuc_dict = {}

    # initialize a count for nucleotides in the sequence
    nuc_count = 0

    # begin looping through VCF file to look for positions that will be replaced
    # and add them to the dictionary
    pos_duplicate_check = 0
    current_position_no_duplicates = 0
    pos_duplicate_count = 0

    # initialize reference name variable. this may not be necessary later
    ref_str = ''

    proc_seq = process_alignment(args.align_file, ref_sequence_list, nuc_count, ref_str)

    last_line = subprocess.check_output(['tail', '-1', args.vcf_file]).decode('UTF-8')
    split_last_line = last_line.split()
    vcf_length = split_last_line[1]
    print(vcf_length)

    vcf_to_dict = process_vcf(args.vcf_file, replace_nuc_dict)
    # print(vcf_to_dict)

    # check_contiguity(vcf_to_dict, vcf_length)

    fixed_unchecked_length_seq = compare_and_fix(vcf_to_dict, ref_sequence_list, output_list, ref_str, nuc_set)

    # assert len(ref_sequence_list) == len(fixed_unchecked_length_seq)
    print("Sequence length after processing: ", len(ref_sequence_list))

    length_fix = ensure_seq_length(ref_sequence_list, vcf_length)

    #Write output sequence
    fasta_output_file = open(args.out_file, 'w')
    fasta_output_file.write(proc_seq)
    fasta_output_file.write("\n")
    fasta_output_file.write(length_fix)
    fasta_output_file.close()




    # # open sequence file that will have the nucleotides replaced
    # seq_file = open(args.align_file, 'r')
    #
    # # read sequence file, split the file on new line characters
    # # and take the first chunk as the taxon name
    # read_seq = seq_file.read()
    # split_seqs = read_seq.split('\n')
    # # tax_name = split_seqs[0:1]
    # tax_name = split_seqs[0]
    # print(tax_name)
    # str_tax_name = ''.join(tax_name)
    #
    # # seperate out the original sequence produced by mpileup
    # seq = split_seqs[1:2]
    # seq_str = seq[0]
    # for nuc in seq_str:
    #     ref_sequence_list.append(nuc)
    #     nuc_count+=1


    # # initialize reference name variable. this may not be necessary later
    # ref_str = ''

    # # initialize dictionary for location of nucleotide to be replaced
    # # and the nucleotide that will go there
    # replace_nuc_dict = {}
    #
    # # begin looping through VCF file to look for positions that will be replaced
    # # and add them to the dictionary
    # pos_duplicate_check = 0
    # current_position_no_duplicates = 0
    # pos_duplicate_count = 0

    # last_line = subprocess.check_output(['tail', '-1', args.vcf_file]).decode('UTF-8')
    # split_last_line = last_line.split()
    # vcf_length = split_last_line[1]
    # print(vcf_length)


    # with open(args.vcf_file) as vcf:
    #     for line_num, line in enumerate(vcf):
    #         if not line.startswith("#"):
    #             splitter = line.split()
    #             # ref_name = splitter[0:1]
    #             ref_name = splitter[0]
    #             # print(ref_name)
    #             ref_str = ''.join(ref_name)
    #             # pos = splitter[1:2]
    #             pos = splitter[1]
    #             # print(pos)
                # if pos_duplicate_check != pos:
                #     pos_duplicate_check = pos
                #
                # elif pos_duplicate_check == pos:
                #     pos_duplicate_count+=1
                #
                # ref = splitter[3:4]
                # alt = splitter[4:5]
                # if 'N' in ref:
                #     if alt in nuc_set:
                #         pos_str = ''.join(pos)
                #         str_alt = ''.join(alt)
                #         split_alt = str_alt.split(',')
                #         alt_nuc = split_alt[0]
                #         replace_nuc_dict[pos_str] = alt_nuc
                #     elif alt not in nuc_set:
                #         pos_str = ''.join(pos)
                #         str_alt = ''.join(alt)
                #         split_alt = str_alt.split(',')
                #         alt_nuc = 'N'
                #         replace_nuc_dict[pos_str] = alt_nuc
                #     else:
                #         print('YOU HAVE A LARGE PROBLEM WITH YOUR VCFFIXER.PY')





    # loop through the dictionary and find the positions in the list
    # that matches the position for the new nucleotide in the dictionary
    # replace the nucleotide in the list at that position
    # with the nucleotide from the dictionary
    # print(len(replace_nuc_dict))
    # print(replace_nuc_dict)

    # print(replace_nuc_dict[vcf_length])
    # find_len_diff = find_length_difference(replace_nuc_dict, vcf_length)
    # print(find_len_diff)
    # find_first_position(replace_nuc_dict, vcf_length)
    # for key, value in replace_nuc_dict.items():
    #     ref_sequence_list[int(key) - 1] = value
    #
    # # turn the list into a string
    # fixed_seq = ''.join(ref_sequence_list)
    #
    # # write the taxon name, a new line character and the fixed sequence to file
    # new_file = open(args.out_file,'w')
    # new_file.write(str_tax_name)
    # new_file.write('\n')
    # new_file.write(fixed_seq)

def check_contiguity(vcf_dict, seq_len):
    """Makes sure all lines are present in the VCF"""
    count_check = 0
    for num in range(0, int(seq_len)):
        count_check+=1
        # print(num)
        check_key = vcf_dict[num + 1]
        assert len(check_key) > 0




def process_vcf(vcf_file, dict):
    """Read VCF file and add nucleotides and positions to a dictionary. \
    Dict protects against multiple nucleotide options at a single position"""
    with open(vcf_file) as vcf:
        for line_num, line in enumerate(vcf):
            if not line.startswith("#"):
                splitter = line.split()
                nucs_list = []
                ref_name = splitter[0]
                # ref_str = ''.join(ref_name)
                pos = splitter[1]
                # ref = splitter[3:4]
                # alt = splitter[4:5]
                ref = splitter[3]
                alt = splitter[4]
                # print(ref)
                # print(alt)

                dupe_check = check_dict_duplicate(dict, pos)
                if dupe_check == False:
                    nucs_list.append(ref)
                    nucs_list.append(alt)
                    # print(nucs_list)
                    dict[int(pos)] = nucs_list

    return dict


def compare_and_fix(vcf_dict, align_list, output_list, taxon_name, nuc_set_):
    """Reads over alignment sequence and checks if an alternative nucleotide was recorded where mpileup placed an N. \
    The N is replaced by the alternative nucleotide"""
    print("align list length prior to adjustment: ", len(align_list))
    # print(align_list)
    for num, nuc in enumerate(align_list):
        # print(nuc)
        if nuc.upper() == "N":
            pos_in_vcf = num + 1
            # ref_and_alt = vcf_dict[pos_in_vcf]
            ref_and_alt = vcf_dict.get(pos_in_vcf)
            if ref_and_alt == None:
                output_list.append("N")
            elif ref_and_alt != None:
                ref = ref_and_alt[0]
                alt = ref_and_alt[1]
                split_alts = alt.split(',')
                if len(split_alts) > 1:
                    # print(split_alts)

                    # OK, here we just take the first alternative nucleotide
                    # There is probably a better way to handle this but for now, this works
                    selected_alt = split_alts[0]
                    if selected_alt in nuc_set_:
                        # output_list.append(selected_alt.upper())
                        align_list[num] = selected_alt.upper()
                        print("N replaced with {} and position {} from vcf position {}.".format(selected_alt, num, pos_in_vcf))

        # elif nuc.upper() in nuc_set_:
        #     output_list.append(nuc.upper())

        elif nuc.upper() != "N":
            if nuc.upper() not in nuc_set_:
            # output_list.append("N")
                align_list[num] = "N"



    #return output_list

def ensure_seq_length(processed_seq, vcf_length_):
    """adds Ns at the end of the sequence to make the new, corrected sequence \
    as long as the length the VCF says it should be. Convert list to string \
    and return to be written to file."""
    print(len(processed_seq))
    print(vcf_length_)
    if len(processed_seq) < int(vcf_length_):
        diff_in_len = int(vcf_length_) - len(processed_seq)
        print(diff_in_len)
        for num in range(0, diff_in_len):
            processed_seq.append("N")
    print(len(processed_seq))
    processed_seq = ''.join(processed_seq)

    return processed_seq



def check_dict_duplicate(dict, check_key):
    """Check if a key/value is already present in a dictionary."""
    # print(dict.keys())
    if check_key in dict.keys():
        print("DUPLICATE POSITION FOUND at position: ", check_key)
        return True
    else:
        return False

def process_alignment(align_file, output_dict, base_count, taxon_name):
    """read alignment file, add individual nucleotides to a list, count the bases and store the taxon name"""
    print("Processing alignment.")
    # open sequence file that will have the nucleotides replaced
    seq_file = open(align_file, 'r')

    # read sequence file, split the file on new line characters
    # and take the first chunk as the taxon name
    read_seq = seq_file.read()
    split_seqs = read_seq.split('\n')
    # tax_name = split_seqs[0:1]
    tax_name = split_seqs[0]
    print(tax_name)
    str_tax_name = ''.join(tax_name)

    # seperate out the original sequence produced by mpileup
    seq = split_seqs[1]
    # seq_str = seq[0]
    # print(seq_str)
    for nuc in seq:
        output_dict.append(nuc)
        base_count+=1

    print("Alignment nucleotide count: ", base_count)

    return tax_name


def find_length_difference(dict_of_positions, vcf_length):
    int_len = int(vcf_length)
    str_vcf_length = str(int_len)
    if str_vcf_length not in dict_of_positions.keys():
        find_length_difference(dict_of_positions, (int_len - 1))
    elif str_vcf_length in dict_of_positions.keys():
        return str_vcf_length

def find_first_position(dict_of_positions, end_len):
    starting_pos = 0
    while starting_pos < int(end_len) and starting_pos < 800:
        if str(starting_pos) in dict_of_positions.keys():
            print("--")
        else:
            starting_pos+=1
            print("NOT FOUND: ", starting_pos)


if __name__ == '__main__':
    main()
