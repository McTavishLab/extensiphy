#! /usr/bin/python

import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_file')
    parser.add_argument('--out_file')
    parser.add_argument('--seq_name')
    return parser.parse_args()


def main():
    args = parse_args()
    
    nucleotides = {"A", "C", "G", "T"}

    accepted_degenerate_bases = set("N") # this needs to be expanded
    new_seq = ''
    pos_tracker = 0
    indel_tracker = False
    with open(args.vcf_file) as vcf:
        for line_num, line in enumerate(vcf):

            if not line.startswith("#"):

                splitter = line.split()
                ref_name = splitter[0:1]
                pos = splitter[1:2]
                ref = splitter[3:4]
                alt = splitter[4:5]
                #print(pos)
                #print(type(ref[0]))
                #print(alt[0])
                #print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                #if int(pos[0]) != pos_tracker + 1:
                #    new_seq = new_seq + ("N" * (int(pos[0]) - 1))

                # Handles adding N's to the sequence if there are missing nucleotides from the start of the vcf
                #if len(new_seq) == 0 and pos > 1:
                #    new_seq = new_seq + ('N' * (int(pos[0]) - 1))

                # handles adding N's to the sequence when a gap is detected in the vcf
                if int(pos[0]) != pos_tracker + 1:
                    new_seq = new_seq + ("N" * ((int(pos[0]) - pos_tracker) - 1))
                    #print(pos)
                    #print(pos_tracker)
                    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~")
                if int(pos[0]) == pos_tracker:
                    print("INDEL FOUND AT POS ", pos[0])
                    indel_tracker = True

                if indel_tracker == False and str(ref[0][0]) == "N":
                    new_seq = new_seq + alt[0][0]
                    print("N FOUND AT THIS POSITION")
                    print(pos[0])
                    print(alt[0])
                    print("~~~~~~~~~")
               
                elif indel_tracker == False and alt[0][0] not in nucleotides:
                    new_seq = new_seq + ref[0][0]

                elif indel_tracker == False and alt[0][0] in nucleotides:
                    new_seq = new_seq + alt[0][0]
                
                elif indel_tracker == True:
                    print("INDEL FOUND, SKIPPING")

                # Handles adding N's if there are gaps in the vcf sequence
                pos_tracker = int(pos[0])
                indel_tracker = False
    
    
    output_file = open(args.out_file, "w")
    output_file.write(">")
    output_file.write(args.seq_name)
    output_file.write("\n")
    output_file.write(new_seq)
    output_file.close()




if __name__ == '__main__':
    main()
