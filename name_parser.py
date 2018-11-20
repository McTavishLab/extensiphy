#! /usr/bin/python

import sys
import re
import argparse
import json
import os
import collections
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description = 'Rename sequences in alignments and trees. Use previous dictionaries of names if available')
    parser.add_argument('-o', action='store_true', help='Use to rename taxa in an alignment and tree based on that alignment. Use with arguments --alignment_file [ALIGNMENT FILE] and --tree_file [TREE_FILE]')
    parser.add_argument('-d', action='store_true', help='Use this to rename read files before input into any program, producing a json dictionary of names')
    parser.add_argument('-u', action='store_true', help='Use with aruguments --readset_dir [DIRECTORY] and --dict_file [FILE] to rename taxa being added to tree. ONLY use if the tree and alignment file have already been renamed and a dictionary of names has been produced')
    parser.add_argument('--tree_file')
    parser.add_argument('--alignment_file')
    parser.add_argument('--dict_file')
    parser.add_argument('--newtaxa_dir')
    return parser.parse_args()


def main():
    # read in files
    args = parse_args()

    if args.d == True:

        # This is EXTREMELY overengineered. So many counts. I'm sorry if you have to edit this
        # ill try to fix it to be an actual sorting algorithm soon
        # i needed results quick! *sob*
        files = os.listdir(args.newtaxa_dir)
        match1_dict = {}
        match2_dict = {}
        sd = sorted(files)
        pair_1_count = 0
        pair_2_count = 0
        file_count = 0
        read_set_count = 0
        for file in (sd):
            file_count+=1
            read_set_count+=1
            if file_count == 1:
                pair_1_count+=1
                match1_dict["taxon_" + str(pair_1_count)] = file
            elif file_count == 2:
                pair_2_count+=1
                match2_dict["taxon_" + str(pair_2_count)] = file
                file_count = 0
        # print(match1_dict)
        # print(match2_dict)

        od1 = collections.OrderedDict(sorted(match1_dict.items()))
        od2 = collections.OrderedDict(sorted(match2_dict.items()))

        new_names_1 = {}
        for key, value in od1.items():
            new_names_1[key + '_R1.fastq'] = value
        print(new_names_1)

        new_names_2 = {}
        for key, value in od2.items():
            new_names_2[key + '_R2.fastq'] = value
        print(new_names_2)

        os.chdir(args.newtaxa_dir)

        for key, value in new_names_1.items():
            os.rename(str(value), str(key))

        for key, value in new_names_2.items():
            os.rename(str(value), str(key))

        with open('taxon_names_dict_set_1.txt', 'w') as f:
            json.dump(new_names_1, f)

        with open('taxon_names_dict_set_2.txt', 'w') as f:
            json.dump(new_names_2, f)







    # reads in a directory of new taxa fastas produced by multi_map.sh
    # and a dictionary of names produced by an earlier run of name_parser.py
    # produces:
    # new seperate fasta files with renamed sequences
    # new updated dictionary for all taxa with new and old names

    elif args.u == True:
        name_dict = {}
        current_OTU_count = 0
        OTU_nums = []
        taxa_dir = os.listdir(args.newtaxa_dir)
        #new_taxa_dict =
        with open(args.dict_file) as json_data:
            read_dict = json_data.read()
            name_dict = json.loads(read_dict)
            for key, value in name_dict.iteritems():
                key = str(key)
                split_key = key.split('_')
                for item in split_key[1:]:
                    OTU_nums.append(int(item))
        current_OTU_count = max(OTU_nums)
        for file in taxa_dir:
            #new_name_file = open('OTU_' + str(current_OTU_count) + '.fas')
            current_OTU_count+=1
            new_name_file = open(args.newtaxa_dir + '/' + 'OTU_' + str(current_OTU_count) + '.fas', 'w')
            look = open(args.newtaxa_dir + '/' + file,'r')
            examine = look.read()
            split_read = examine.split('\n')
            taxon_name = split_read[0]
            name_dict['OTU_' + str(current_OTU_count)] = taxon_name
            for key, value in name_dict.iteritems():
                if re.match(value, examine):
                    replaces = re.sub(value, key, examine)
                    new_name_file.write('>')
                    new_name_file.write(replaces)
                    new_name_file.write('\n')
        OTU_name_dict = "OTU_name_dict.txt"
        # write dictionary of names file
        with open(args.dict_file, 'w') as f:
            json.dump(name_dict, f)


    # takes in:
    # alignment and tree tree files
    # produces:
    # renamed alignment file, renamed tree file and dictionary of old and new names
    elif args.o == True:
        alignment = open(args.alignment_file)
        tree = open(args.tree_file)

        readdata = alignment.read()

        read_topo = tree.read()

        # split on newlines for alignment data
        strp = readdata.strip()
        splt = strp.split('>')

        # drop wonky space at beginning of list of splits
        clean_splt = splt[1:]

        # initiate dictionary to contain old and new names
        newnames_oldnames = {}
        OTU_count = 0

        # iterate over sequences and add original names plus new names to dictionary
        # keys = new newnames
        # values = old names
        for i in clean_splt:
            count = 0
            OTU_count+=1
            j = i.split('\n')
            for segment in j:
                count+=1
                if count==1:
                    newnames_oldnames['OTU_' + str(OTU_count)] = segment

        # open file to output dictionary
        OTU_name_dict = "OTU_name_dict.txt"
        new_alignment = open("renamed_alignment.fa", 'w')

        # iterate through dictionary and list of sequences
        # replace old names with new names
        for key, value in newnames_oldnames.iteritems():
            for item in clean_splt:
                if re.match(value, item):
                    replaces = re.sub(value, key, item)
                    new_alignment.write(">")
                    new_alignment.write(replaces)
                    new_alignment.write('\n')

        # read tree file
        # topology = open(tree, 'r')
        # read_topo = topology.read()
        tree_name_count = 0

        # initialize str to be used as tree saving device outside of for loop
        swap = ''
        swap = re.sub(value, key, read_topo)

        # iterate over tree and replace old names with new names
        for key, value in newnames_oldnames.iteritems():
            swap = re.sub(value, key, swap)

        # write new tree file
        new_tree = open('renamed_tree.tre','w')
        new_tree.write(swap)

        # write dictionary of names file
        with open(OTU_name_dict, 'w') as f:
            json.dump(newnames_oldnames, f)

            #for key, value in newnames_oldnames.items():
                #f.write('%s:%s\n' % (key, value))


        new_tree.close()
        new_alignment.close()

if __name__ == '__main__':
    main()
