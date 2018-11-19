#! /usr/bin/python

import sys
import re
import argparse
import json
import os

def parse_args():
    parser = argparse.ArgumentParser(description = 'Rename sequences in alignments and trees. Use previous dictionaries of names if available')
    parser.add_argument('--rename_existing', action='store_true', help='Use to rename taxa in an alignment and tree based on that alignment')
    parser.add_argument('--tree_file')
    parser.add_argument('--alignment_file')
    parser.add_argument('--provide_names', action='store_true', help='Use this option in conjuction with --dict_file to provide a dictionary of names from a previous phycorder run')
    parser.add_argument('--dict_file')
    parser.add_argument('--newtaxa_dir')
    parser.add_argument('--rename_taxon', action='store_true', help='Use with aruguments --readset_dir to rename taxa being added to tree')
    return parser.parse_args()


def main():
    # read in files
    args = parse_args()

    if args.provide_names == True and args.rename_taxon == False:
        alignment = open(args.alignment_file)
        tree = open(args.tree_file)

        readdata = alignment.read()
        read_topo = tree.read()

        with open(args.dict_file) as json_data:
            read_dict = json_data.read()
            name_dict = json.loads(read_dict)
            for key, value in name_dict.iteritems():
                key = str(key)
                value = str(value)
                read_topo = str(read_topo)
                if key in read_topo:
                #if re.match(key, read_topo):
                    print("found a correct OTU name in tree file")
                elif value in read_topo:
                #elif re.match(value, read_topo):
                     print("found old name in tree file")



    elif args.rename_taxon == True and args.provide_names == True:
        name_dict = {}
        current_OTU_count = 0
        OTU_nums = []
        taxa_dir = os.listdir(args.newtaxa_dir)
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
            look = open(args.newtaxa_dir + '/' + file,'r')
            examine = look.read()
            






        # # split on newlines for alignment data
        # strp = readdata.strip()
        # splt = strp.split('>')
        #
        # # drop wonky space at beginning of list of splits
        # clean_splt = splt[1:]
        #
        # # initiate dictionary to contain old and new names
        # newnames_oldnames = {}
        # OTU_count = 0
        #
        # # iterate over sequences and add original names plus new names to dictionary
        # # keys = new newnames
        # # values = old names
        # for i in clean_splt:
        #     count = 0
        #     OTU_count+=1
        #     j = i.split('\n')
        #     for segment in j:
        #         count+=1
        #         if count==1:
        #             newnames_oldnames['OTU_' + str(OTU_count)] = segment
        #
        # # open file to output dictionary
        # OTU_name_dict = "OTU_name_dict.txt"
        # new_alignment = open("renamed_alignment.fa", 'w')
        #
        # # iterate through dictionary and list of sequences
        # # replace old names with new names
        # for key, value in newnames_oldnames.iteritems():
        #     for item in clean_splt:
        #         if re.match(value, item):
        #             replaces = re.sub(value, key, item)
        #             new_alignment.write(">")
        #             new_alignment.write(replaces)
        #             new_alignment.write('\n')
        #
        # # read tree file
        # # topology = open(tree, 'r')
        # # read_topo = topology.read()
        # tree_name_count = 0
        #
        # # initialize str to be used as tree saving device outside of for loop
        # swap = ''
        # swap = re.sub(value, key, read_topo)
        #
        # # iterate over tree and replace old names with new names
        # for key, value in newnames_oldnames.iteritems():
        #     swap = re.sub(value, key, swap)
        #
        # # write new tree file
        # new_tree = open('renamed_tree.tre','w')
        # new_tree.write(swap)
        #
        # # write dictionary of names file
        # with open(OTU_name_dict, 'w') as f:
        #     for key, value in newnames_oldnames.items():
        #         f.write('%s:%s\n' % (key, value))
        #
        #
        # new_tree.close()
        # new_alignment.close()

    else:
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
