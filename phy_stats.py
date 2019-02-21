#! /usr/bin/python3

import sys
import os
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_file')
    parser.add_argument('--taxon_output_file')
    parser.add_argument('--align_output_file')
    return parser.parse_args()

def main():
    args = parse_args()

    degen = ['I', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D']
    nucleotides = ['A', 'C', 'G', 'T']
    gaps = ['-']
    Ns = ['N']

    align_n_count=0
    align_degen_count=0
    align_gap_count=0
    align_nucleotide_count=0
    taxon_info_dict = {}

    file_segments = []
    with open(args.align_file, 'r') as fasta:

        read_fasta = fasta.read()
        split_file = read_fasta.split('>')
        useful_align = split_file[1:]
        for taxon in useful_align:

            str_aln = ''.join(taxon)
            taxon_name = str_aln.split('\n', 1)[0]
            seq = str_aln.split('\n', 1)[1]
            one_line_seq = "".join(seq.splitlines())
            taxon_info_dict[taxon_name] = {}

            for species in taxon_info_dict:
                if taxon_name == species:
                    taxon_nucleotide_count = 0
                    taxon_gap_count = 0
                    taxon_degen_count = 0
                    taxon_n_count = 0

                    for nuc in one_line_seq:
                        nuc = nuc.upper()
                        #print(nuc)
                        if nuc in nucleotides:
                            align_nucleotide_count+=1
                            taxon_nucleotide_count+=1
                        elif nuc in gaps:
                            align_gap_count+=1
                            taxon_gap_count+=1
                        elif nuc in degen:
                            align_degen_count+=1
                            taxon_degen_count+=1
                        elif nuc in Ns:
                            align_n_count+=1
                            taxon_n_count+=1



            taxon_info_dict[species]['taxon nucleotide count'] = taxon_nucleotide_count
            taxon_info_dict[species]['taxon gap count'] = taxon_gap_count
            taxon_info_dict[species]['taxon degeneracy count'] = taxon_degen_count
            taxon_info_dict[species]['taxon N count'] = taxon_n_count

    # taxon_info_dict['Total Alignment Info'] = {}
    # taxon_info_dict['Total Alignment Info']['alignment nucleotide count'] = align_nucleotide_count
    # taxon_info_dict['Total Alignment Info']['alignment gap count'] = align_gap_count
    # taxon_info_dict['Total Alignment Info']['alignment degeneracy count'] = align_degen_count
    # taxon_info_dict['Total Alignment Info']['alignment N count'] = align_n_count

    alignment_info_dict = {}
    alignment_info_dict['Alignment nucleotide count'] = align_nucleotide_count
    alignment_info_dict['Alignment gap count'] = align_gap_count
    alignment_info_dict['Alignment degeneracy count'] = align_degen_count
    alignment_info_dict['Alignment N count'] = align_n_count

    # convert nested dictionary into a pandas dataframe
    tax_info_df = pd.DataFrame.from_dict(taxon_info_dict, orient = 'index')
    aln_info_df = pd.DataFrame.from_dict(alignment_info_dict, orient = 'index')

    tax_info_df.to_csv(args.taxon_output_file, sep = '\t' )
    aln_info_df.to_csv(args.align_output_file, sep = '\t' )

if __name__ == '__main__':
    main()
