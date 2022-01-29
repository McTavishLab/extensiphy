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

    # create an empty list that will hold the fixed sequence

    cns_list, tax_name = process_alignment(args.align_file)
    cns_seq_len = len(cns_list)
    # print(cns_list)
    print(cns_seq_len)

    last_line = subprocess.check_output(['tail', '-1', args.vcf_file]).decode('UTF-8')
    last_vcf_pos = last_line.split()[1]
    last_vcf_base = last_line.split()[2]

    if int(last_vcf_pos) > int(cns_seq_len):
        assert(last_vcf_base=='N')
        print("length difference detected")
        length_fix = ensure_seq_length(cns_list, last_vcf_pos)
    else:
        length_fix = ''


    vcf_dict = process_vcf(args.vcf_file)
    # print(vcf_dict)


    # check_contiguity(vcf_to_dict, vcf_length)
    if _DEBUG:
        print("VCFFIXER process VCF complete <><><>")

    output_list = compare_and_fix(vcf_dict, cns_list)

    # assert len(ref_sequence_list) == len(fixed_unchecked_length_seq)
    if _DEBUG:
        print("VCFFIXER Sequence length after processing: ", len(output_list))
        print("VCFFIXER Writing to file. <><><>")

    #Write output sequence
    fasta_output_file = open(args.out_file, 'w')
    fasta_output_file.write(tax_name+"\n")
    for i, pos in enumerate(output_list):
        i+=1
        fasta_output_file.write(pos)
        # if i%140 == 0:
        #     fasta_output_file.write('\n')
    fasta_output_file.write(length_fix)
    fasta_output_file.close()


def check_contiguity(vcf_dict, seq_len):
    """Makes sure all lines are present in the VCF"""
    count_check = 0
    for num in range(0, int(seq_len)):
        count_check+=1
        # print(num)
        check_key = vcf_dict[num + 1]
        assert len(check_key) > 0




def process_vcf(vcf_file):
    """Read VCF file and add nucleotides and positions to a dictionary. \
    Dict protects against multiple nucleotide options at a single position"""
    vcf_dict = {}
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
                if pos in vcf_dict.keys():
                    pass
                    #what should happen if position is a duplicate?
                else:
                    nucs_list.append(ref)
                    nucs_list.append(alt)
                    # print(nucs_list)
                    vcf_dict[int(pos)] = nucs_list
    return vcf_dict


def compare_and_fix(vcf_dict, cns_list):
    """Reads over alignment sequence and checks if an alternative nucleotide was recorded where mpileup placed an N. \
    The N is replaced by the alternative nucleotide"""
    if _DEBUG:
        print("align list length prior to adjustment: ", len(cns_list))
    output_list = []
    for num, nuc in enumerate(cns_list):
        # print(nuc)
        if nuc.upper() == "N":
            pos_in_vcf = num + 1
            ref_and_alt = vcf_dict.get(pos_in_vcf)
            if ref_and_alt == None:
                output_list.append(nuc)
            elif ref_and_alt != None:
                ref = ref_and_alt[0]
                if ref != 'N':
                    output_list.append(nuc)
                else:
                    alt = ref_and_alt[1]
                    split_alts = alt.split(',')
                    if len(split_alts) > 1:
                        # print(split_alts)
                        # OK, here we just take the first alternative nucleotide
                        # There is probably a better way to handle this but for now, this works
                        selected_alt = split_alts[0]
                    else:
                        selected_alt = alt
                    if selected_alt in ['A', 'C', 'G', 'T']:
                        output_list.append(selected_alt.upper())
                        if _DEBUG:
                            print("N replaced with {} and position {} from vcf position {}.".format(selected_alt, num, pos_in_vcf))
                    else:
                        output_list.append(nuc)
                        ## Q. What if it's not one of those bases??
        else:
            output_list.append(nuc)
    assert(len(output_list) == len(cns_list)), "{},{}".format(len(output_list), len(cns_list))
    return output_list


def ensure_seq_length(cns_dict, last_vcf_pos):
    """adds Ns at the end of the sequence to make the new, corrected sequence \
    as long as the length the VCF says it should be. Convert list to string \
    and return to be written to file."""
    if _DEBUG:
        print(len(processed_seq))
        print(vcf_length)
    if len(processed_seq) < int(last_vcf_pos):
        diff_in_len = int(last_vcf_pos) - len(processed_seq)
        if _DEBUG:
            print(diff_in_len)
        for num in range(0, diff_in_len):
            processed_seq.append("N")
    print(len(processed_seq))
    processed_seq = ''.join(processed_seq)
    return processed_seq


def process_alignment(align_file):
    """read consensus sequence file, add individual nucleotides to a dictionary"""
    if _DEBUG:
        print("Processing alignment.")
    # open sequence file that will have the nucleotides replaced
    seq_file = open(align_file, 'r')

    # read sequence file, split the file on new line characters
    # and take the first chunk as the taxon name
    read_seq = seq_file.read()
    split_seqs = read_seq.split('\n')
    # print(split_seqs)
    # tax_name = split_seqs[0:1]
    tax_name = split_seqs[0]
    if _DEBUG:
        print(tax_name)
    # seperate out the original sequence produced by mpileup
    seq = split_seqs[1]
#    assert((set(seq.upper()))==set('ATGCN')), set(seq.upper())
    # seq_str = seq[0]
    # print(seq_str)
    cns_list = list(seq)
    if _DEBUG:
        print("Alignment nucleotide count: ", len(cns_list))
    return cns_list, tax_name


if __name__ == '__main__':
    main()
