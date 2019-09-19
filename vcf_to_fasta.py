#! /usr/bin/python

import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_file')
    parser.add_argument('--out_file')
    return parser.parse_args()


def main():
    args = parse_args()
    
    nucleotides = {"A", "C", "G", "T"}

    accepted_degenerate_bases = set("N") # this needs to be expanded
    new_seq = ''
    pos_tracker = 0
    with open(args.vcf_file) as vcf:
        for line_num, line in enumerate(vcf):

            if not line.startswith("#"):

                splitter = line.split()
                ref_name = splitter[0:1]
                pos = splitter[1:2]
                ref = splitter[3:4]
                alt = splitter[4:5]
                
                #if int(pos[0]) != pos_tracker + 1:
                #    new_seq = new_seq + ("N" * (int(pos[0]) - 1))

                # Handles adding N's to the sequence if there are missing nucleotides from the start of the vcf
                if len(new_seq) == 0 and pos > 1:
                    new_seq = new_seq + ('N' * (int(pos[0]) - 1))

                # handles adding N's to the sequence when a gap is detected in the vcf
                if int(pos[0]) != pos_tracker + 1:
                    new_seq = new_seq + ("N" * (int(pos[0]) - pos_tracker))
                
                # Handles adding the correct ref or alt nucleotide to the sequence
                if alt[0] not in nucleotides:
                    new_seq = new_seq + ref[0]

                elif alt[0] in nucleotides:
                    new_seq = new_seq + alt[0]
                
                # Handles adding N's if there are gaps in the vcf sequence
                pos_tracker = int(pos[0])

    print(new_seq)
                #print(alt)
                #print("~~~~~~~~~~~~~~~")
    




if __name__ == '__main__':
    main()
