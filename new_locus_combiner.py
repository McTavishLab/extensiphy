#! /usr/bin/python

import os
import argparse
import re
import json
from collections import defaultdict
import csv

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--msa_folder')
    parser.add_argument('--out_file')
    parser.add_argument('--position_csv_file')
    parser.add_argument('--suffix')
    parser.add_argument('--len_filter')
    return parser.parse_args()


def main():
    args = parse_args()
    len_filt = int(args.len_filter)
    assert(type(len_filt) == int)
    dir_of_aligns = args.msa_folder
    # check if path ends with a "/"
    if dir_of_aligns.endswith("/"):
        print("Pathing acceptable")

    else:
        print("Fixing pathing")
        dir_of_aligns = dir_of_aligns + "/"
        print(dir_of_aligns)
    print(dir_of_aligns)
    #msa_list = os.listdir(args.msa_folder)
    msa_list = os.listdir(dir_of_aligns)
    print(msa_list)
   
    name_list = []
    file_name_and_seq_len_dict = {}
    taxon_name_and_seqs = defaultdict(list)
    final_dict = {}
    tuple_name_and_seq_len_list = []
    loci_count = 0
    csv_header = ['locus_position_number', 'locus_file_name', 'locus_length']
    tuple_name_and_seq_len_list.append(csv_header)

    for file_name in msa_list:
        if file_name.endswith(args.suffix):
            locus_name = file_name
            locus_len = 0
            loci_count+=1
            open_file = open(dir_of_aligns + file_name, "r")
            read_file = open_file.read()
            split_file = read_file.split(">")
            for name_and_seq in split_file:
                if len(name_and_seq) > 1:
                    name_seq_split = name_and_seq.split("\n", 1)
                    name = name_seq_split[0]
                    seq = name_seq_split[1]
                    seq = seq.replace("\n","")
                    seq_len = len(seq)
                    print(seq_len)
                    if seq_len >= len_filt:
                        taxon_name_and_seqs[name].append(seq)
                        file_name_and_seq_len_dict[file_name] = seq_len
                        locus_len = seq_len
                        #tuple_name_and_seq_len_list.append(name_and_len)
                        if file_name not in name_list:
                            name_list.append(file_name)
                    elif seq_len < len_filt:
                        print("locus not passing length filter!")
                    #file_name_and_seq_len_dict[file_name] = seq_len
                    #locus_len = seq_len
                    #tuple_name_and_seq_len_list.append(name_and_len)
                    #if file_name not in name_list:
                    #    name_list.append(file_name)
            if seq_len >= len_filt:
                name_and_len = [loci_count, file_name, seq_len]
                tuple_name_and_seq_len_list.append(name_and_len)
    print(tuple_name_and_seq_len_list)
    #print(file_name_and_seq_len_dict)
    #print(taxon_name_and_seqs)
    
    myFile = open(args.position_csv_file, 'w')
    with myFile:
        writer = csv.writer(myFile)
        writer.writerows(tuple_name_and_seq_len_list)
     
    print("Writing complete") 
    
    
    
    
    for num, f_name in enumerate(name_list):
        if f_name in file_name_and_seq_len_dict.keys():
            temp_dict = {}
            temp_dict[file_name_and_seq_len_dict[f_name]] = f_name
            final_dict[num] = temp_dict
            #final_dict[num] = "{" + str(file_name_and_seq_len_dict[f_name]) + ":" + f_name + "}"

    #print("FINAL DICT OF POSITIONS")
    #print(final_dict)

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
        
                        
    





    #with open(args.position_dict_file, 'w') as dict_output:
    #    json.dump(final_dict, dict_output)
#    concat_file.close()



if __name__ == '__main__':
    main()
