#! /usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import dendropy

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree_file')
    return parser.parse_args()

def split_organiser(splits_with_len_dict):
    #COUNT TAXA IN DATASET
    taxa_count = 0
    for split, length in splits_with_len_dict.items():
        for position in split:
            taxa_count+=1

    split_sort_dict = {}

    for split, length in splits_with_len_dict.items():
        taxa_in_split = 0
        for position in split:
            if position == '1':
                taxa_in_split+=1
        if taxa_in_split not in split_sort_dict:
            split_sort_dict[taxa_in_split] = []
            split_sort_dict[taxa_in_split].append(split)
        
        else:
            split_sort_dict[taxa_in_split].append(split)
    
    return split_sort_dict

#CONSTRUCT INDEX OF KEYS IN ORGANIZED SPLIT_SORT_DICT AND SORT BY DECENDING SIZE
def split_indexer(organized_split_dict):
    index = []
    for key, value in organized_split_dict.items():
        index.append(key)
    sorted_index = index[::-1]
    return sorted_index

def nested_split_constructor(index, organized_splits):
    split_count_tracker = 1
    for num in index:
        num_of_splits = len(organized_splits[num])
        prev_split = ''
        current_split = ''
        if split_count_tracker == 1:
            prev_split = organized_splits[num][0]
        elif split_count_tracker > 1:
            print(prev_split) 

        split_count_tracker+=1
            


def main():
    args = parse_args()
    taxa = dendropy.TaxonNamespace()
    mle = dendropy.Tree.get(path=args.tree_file, schema='newick', taxon_namespace=taxa)
    mle_len = mle.length()
    mle.encode_bipartitions()
    pdc = mle.phylogenetic_distance_matrix()
    long_branch_cutoff = 0.4
    
    #for i, t1 in enumerate(mle.taxon_namespace[:-1]):
    #    for t2 in mle.taxon_namespace[i+1:]:
    #        print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, pdc(t1, t2)))

    tax_list = [] 
    print(taxa)
    for taxon in taxa:
        tax_list.append(taxon)

    second_splits_branch_length = 0.0
    second_splits_branch_length_percent_of_tree = 0.0
    grouped_splits = []
    split_n_lens = {}
    total_taxa = 0
    split_encode = mle.bipartition_encoding
    for split in split_encode:
        taxa_in_split_count = 0
        split_branch_len = mle.bipartition_edge_map[split].length
        bipart = mle.bipartition_edge_map[split].bipartition
        split_edge_ratio = float(split_branch_len) / float(mle_len)
        second_splits_branch_length+=mle.bipartition_edge_map[split].length
        second_splits_branch_length_percent_of_tree+=split_edge_ratio
        #print(type(bipart))
        str_bipart = str(bipart)
        #print(str_bipart)
        #print(type(str_bipart))
        split_n_lens[str_bipart] = split_edge_ratio
        grouped_splits.append(str_bipart)
        for taxa in str_bipart:
            #print(taxa)
            taxa_in_split_count+=1
        total_taxa = taxa_in_split_count
    #print(split_n_lens)   

    sort_splits = split_organiser(split_n_lens)
    #print(sort_splits)

    index = split_indexer(sort_splits)
    #print(index)

    nest_splits = nested_split_constructor(index, sort_splits)


if __name__ == '__main__':
    main()
