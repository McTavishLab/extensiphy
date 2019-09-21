#! /usr/bin/python

import os
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_file')
    parser.add_argument('--out_file')
    parser.add_argument('--locus_id')
    parser.add_argument('--locus_size')
    return parser.parse_args()


def main():
    args = parse_args()
    
    
    full_sequences = {}
    locus_chunk = ''
    stored_names = []
    stored_ids = []
    locus_id = "(>1:\d+-\d+\s+(\+|\-)\s" + args.locus_id + "\s(.|\n)+?)="
    print(locus_id)
    locus_id_grabber = re.compile(locus_id, re.S)
    file_name = "##SequenceFile(.*?)\n"
    file_name_grabber = re.compile(file_name)
    
    alignment = open(args.align_file, "r")
    read_align = alignment.read()
    #lines = alignment.readlines()
    locus_search = re.search(locus_id_grabber, read_align)
    if locus_search:
        #print(locus_search.group())
        locus_chunk = locus_search.group()
    #print(locus_chunk)
    name_search = re.findall(file_name_grabber, read_align)
    for num, name in enumerate(name_search):
        locus_num = ">" + str(num+1) + ":.*"
        locus_num_grabber = re.compile(locus_num)
        locus_chunk = re.sub(locus_num_grabber, ">" + name , locus_chunk)
    locus_chunk = re.sub("=", "", locus_chunk)
    locus_chunk = re.sub("> ", ">", locus_chunk)
    #locus_file = open(args.out_file, "w")
    #locus_file.write(locus_chunk)
    #locus_file.close()
    
    seq_len_count = 0
    split_taxa = locus_chunk.split(">")
    for item in split_taxa:
        split_seq_and_name = item.split("\n", 1)
        seq_len = len(str(split_seq_and_name[1:2]))
        #for nuc in seq_len:
        #    print(nuc)
        if seq_len >= seq_len_count:
            seq_len_count = seq_len
    print(seq_len_count)
    if int(seq_len_count) >= int(args.locus_size):
        locus_file = open(args.out_file, "w")
        locus_file.write(locus_chunk)
        locus_file.close()

    elif int(seq_len_count) < int(args.locus_size):
        print("locus", args.locus_id, "was not long enough for inclusion.\n")



    alignment.close()




if __name__ == '__main__':
    main()
