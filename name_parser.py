#! /usr/bin/python

import sys
import re
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description = 'Rename sequences in alignments and trees. Use previous dictionaries of names if available')
    parser.add_argument('--tree_file')
    parser.add_argument('--alignment_file')
    parser.add_argument('--dictionary_file')
    return parser.parse_args()


def main():
    # read in files
    args = parse_args()
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
        for key, value in newnames_oldnames.items():
            f.write('%s:%s\n' % (key, value))


    new_tree.close()
    new_alignment.close()

if __name__ == '__main__':
    main()
