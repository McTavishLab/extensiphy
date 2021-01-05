#! /usr/bin/python
# this program will take in a single concatenate sequence produced by Phycorder
# and a positional file that keeps track of the start location of each
# locus in the concatenated Phycorder file
# it will then seperate each locus into its seperate fasta file


import os
import argparse
#import re
#import json
import csv

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out_file_dir')
    parser.add_argument('--position_csv_file')
    parser.add_argument('--concatenated_fasta')
    return parser.parse_args()


def main():
    args = parse_args()

    pos_list = []
    loci_count = 0
    #nuc_count = 0
    #json_data = open(args.position_dict_file,'r')
    #read_dict = json_data.read()
    #pos_dict = json.loads(read_dict)
    with open(args.position_csv_file) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
#            print(row['locus_position_number'], row['locus_file_name'], row['locus_length'])
            pos_list.append(row)
#            loci_count+=1
#            print(loci_count)
    print(pos_list)


#    print(pos_dict)

    #seq_file = open(args.concatenated_fasta, 'r')
    #read_seq = seq_file.read()
    #split_file = read_seq.split('>')

    #print(split_file)

    os.mkdir(args.out_file_dir)

    seq_file = open(args.concatenated_fasta, 'r')
    read_seq = seq_file.read()
    split_file = read_seq.split('>')
    #print(split_file)

    start_pos = 0
    loci_starts = []
    loci_starts.append(start_pos)

#    for seq_number, length_and_loci_name in pos_dict.items():
#        loci_count+=1

    locus_number = 0
    for locus_dict in pos_list:
#        print(locus_dict)
        loci_count+=1
        locus_num_check = int(locus_dict['locus_position_number'])
        print(locus_num_check)
        assert locus_num_check == locus_number + 1
        locus_end_pos = int(locus_dict['locus_length'])
        loci_starts.append(locus_end_pos + start_pos)
        start_pos = start_pos + locus_end_pos
        locus_number = locus_num_check
#        print(start_pos)
    print(loci_starts)

    start_stop_pos_list = []
    print(locus_number)
    count = 0
    for num in loci_starts:
        if count < locus_number:
            locus_name = pos_list[count]['locus_file_name']
            #print(locus_name)
            start = num
            count+=1
            stop = loci_starts[count]
            start_stop_tuple = (start, stop)
            start_stop_pos_list.append(start_stop_tuple)
    print(start_stop_pos_list)

    count_two = 0
    for start_stop_tuple in start_stop_pos_list:
        locus_specific_taxon_and_seq_dict = {}
        #print(locus_specific_taxon_and_seq_dict)
        locus_name = pos_list[count_two]['locus_file_name']
        open_new_file = open(args.out_file_dir + "/" + locus_name,'w')
        for taxon_and_seq in split_file:
            if len(taxon_and_seq) > 1:
                #print(taxon_and_seq)
                split_taxon_and_seq = taxon_and_seq.split("\n", 1)
                taxon_name = split_taxon_and_seq[0]
                taxon_seq = split_taxon_and_seq[1]
                taxon_seq = taxon_seq.replace("\n","")
                #print(taxon_seq)
                seq_chunk = taxon_seq[start_stop_tuple[0]:start_stop_tuple[1]]
                #print(seq_chunk)
                locus_specific_taxon_and_seq_dict[taxon_name] = seq_chunk
        #print(locus_specific_taxon_and_seq_dict)
        count_two+=1
        for tax_name, tax_seq in locus_specific_taxon_and_seq_dict.items():
            #print(tax_seq)
            open_new_file.write(">")
            open_new_file.write(tax_name)
            open_new_file.write("\n")
            open_new_file.write(tax_seq)
            open_new_file.write("\n")
        open_new_file.close()





if __name__ == '__main__':
    main()