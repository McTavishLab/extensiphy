#! /usr/bin/python
# this program will take in a single concatenate sequence produced by Phycorder
# and a positional file that keeps track of the start location of each locus in the concatenated Phycorder file
# it will then seperate each locus into its seperate fasta file


import os
import argparse
import re
import json

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_file_dir')
    parser.add_argument('--position_dict_file')
    parser.add_argument('--concatenated_fasta')
    return parser.parse_args()


def main():
    args = parse_args()
    
    
    loci_count = 0
    nuc_count = 0
    json_data = open(args.position_dict_file,'r')
    read_dict = json_data.read()
    pos_dict = json.loads(read_dict)
    print(pos_dict)

    seq_file = open(args.concatenated_fasta, 'r')
    read_seq = seq_file.read()
    split_file = read_seq.split('\n')
    #print(split_file)

    start_pos = 0
    loci_starts = []
    loci_starts.append(start_pos)
    
    for seq_number, length_and_loci_name in pos_dict.items():
        loci_count+=1

    for i in range(0, loci_count):
        seq_info = pos_dict[str(i)]
        seq_end = list(seq_info.keys())[0]
        print(seq_end)
        start_pos = start_pos + int(seq_end)
        loci_starts.append(start_pos)
    print(loci_starts)


    


    
    #for seq_number, length_and_loci_name  in pos_dict.items():
        
    #    if loci_count == 0:
            
    #        loci_starts.append(start_pos)
    #    for leng, name in length_and_loci_name.items():
            
    #        start_pos = start_pos + int(leng)
    #        loci_starts.append(start_pos)
    #    loci_count+=1
        #print(loci_count)

    #print(loci_starts)

    second_loci_count = 0
    for num in loci_starts:
        second_loci_count+=1
        #print(second_loci_count)
        #print(loci_count)
        if second_loci_count == loci_count:
            seq = split_file[1][num:]
            #print(seq)
            
            seq_info = pos_dict[str(second_loci_count - 1)]
            print(seq_info)
            loci_name = list(seq_info.values())[0]
            print(loci_name)
            output = open(str(args.out_file_dir) + '/' + str(split_file[0]).replace(">","") + str(loci_name), "w")
            output.write(split_file[0])
            output.write("\n")
            output.write(seq)
        
        elif second_loci_count < loci_count:
            seq = split_file[1][num:loci_starts[second_loci_count]]
            #print(seq)
            #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            #print(second_loci_count)
            seq_info = pos_dict[str(second_loci_count - 1)]
            print(seq_info)
            loci_name = list(seq_info.values())[0]
            print(loci_name)
            output = open(str(args.out_file_dir) + "/" + str(split_file[0]).replace(">","") + str(loci_name), "w")
            output.write(split_file[0])
            output.write("\n")
            output.write(seq)

if __name__ == '__main__':
    main() 




