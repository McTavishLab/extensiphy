#! /usr/bin/python

import os
import argparse
import re
import json
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--msa_folder')
    parser.add_argument('--out_file')
    parser.add_argument('--position_dict_file')
    parser.add_argument('--suffix')
    return parser.parse_args()


def main():
    args = parse_args()

    msa_list = os.listdir(args.msa_folder)
    print(msa_list)
   
    name_list = []
    file_name_and_seq_len_dict = {}
    taxon_name_and_seqs = defaultdict(list)
    final_dict = {}

    loci_count = 0
    for file_name in msa_list:
        if file_name.endswith(args.suffix):
            loci_count+=1
            open_file = open(args.msa_folder + file_name, "r")
            read_file = open_file.read()
            split_file = read_file.split(">")
            for name_and_seq in split_file:
                if len(name_and_seq) > 1:
                    name_seq_split = name_and_seq.split("\n", 1)
                    name = name_seq_split[0]
                    seq = name_seq_split[1]
                    seq = seq.replace("\n","")
                    seq_len = len(seq)
                    #print(seq_len)
                    taxon_name_and_seqs[name].append(seq)
                    file_name_and_seq_len_dict[file_name] = seq_len
                    if file_name not in name_list:
                        name_list.append(file_name)


    #print(file_name_and_seq_len_dict)
    #print(taxon_name_and_seqs)
    
    for num, f_name in enumerate(name_list):
        if f_name in file_name_and_seq_len_dict.keys():
            temp_dict = {}
            temp_dict[file_name_and_seq_len_dict[f_name]] = f_name
            final_dict[num] = temp_dict
            #final_dict[num] = "{" + str(file_name_and_seq_len_dict[f_name]) + ":" + f_name + "}"

    print("FINAL DICT OF POSITIONS")
    print(final_dict)

    msa_file = open(args.out_file,'w')
    for name, seqs in taxon_name_and_seqs.items():
        #print(name)
        #print(seqs)
        joined_seqs = ''.join(seqs)
        msa_file.write(">")
        msa_file.write(name)
        msa_file.write("\n")
        msa_file.write(joined_seqs)
        msa_file.write("\n")

    msa_file.close()
        
                        
    





    with open(args.position_dict_file, 'w') as dict_output:
        json.dump(final_dict, dict_output)
#    concat_file.close()



if __name__ == '__main__':
    main()
