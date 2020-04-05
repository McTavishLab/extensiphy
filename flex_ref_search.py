#! /usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import statistics
import numpy
import dendropy

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', help='phylogenetic tree file in newick format')
    parser.add_argument('--dist', default='long', help='branch length distance leading to references (DEFAULT: long)(Options: short, long')
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

#CONSTRUCT INDEX OF KEYS IN ORGANIZED SPLIT_SORT_DICT AND SORT BY SMALLEST TO LARGEST
def split_indexer(organized_split_dict):
    index = []
    for key, value in organized_split_dict.items():
        index.append(key)
    #sorted_index = index[::-1]
    index.sort()
    #return sorted_index
    return index


def combine_splits(split_1, split_2):
    list_of_combined_taxa_positions = list(map(list, zip(split_1, split_2)))

    return list_of_combined_taxa_positions

#MAIN FUNCTION TO ITERATE THROUGH SPLITS AND DECIDE WHAT SPLITS ARE INCLUDED FOR FINDING REFERENCES
def nested_split_constructor(index, organized_splits, splits_and_ratios_dict, cutoff):
    final_splits_list = []
    #print(cutoff)
    num_count = 0
    for num in index:
        num_count+=1
        if num != 0:
            print(num)
            for split in organized_splits[num]:
                #print(splits_and_ratios_dict[split])
                lookat_split = splits_and_ratios_dict[split]
                if lookat_split >= cutoff:
                    #print("found")
                    final_splits_list.append(split)
                
                elif lookat_split < cutoff:
                    if num != max(index):
                        next_split_set = organized_splits[index[num_count]]
                        print(split)
                        print(next_split_set)
                    elif num == max(index):
                       print("waffle") 
    
    
    return final_splits_list


def get_dists(taxa_list, phylo_distance_matrix):
    finish = 0
    dist_list = []
    avgs_list = []
    single_taxa = []
    best_dist = 0
    best_ref = ''
    for subject_tax in taxa_list:
        best_avgs_for_subject = []
        if len(taxa_list) == 1:
            best_ref = subject_tax
            finish = 1
            #print("single best ref")
            return best_ref
        elif len(taxa_list) != 1:
            for selected_tax in taxa_list:
                if selected_tax != subject_tax:
                    dists = phylo_distance_matrix(subject_tax, selected_tax)
                    best_avgs_for_subject.append(dists)
            avg = numpy.mean(best_avgs_for_subject)
            avgs_list.append(avg)
    
    if finish != 1:
        #print("multiple options for best ref, narrowing down")
        best_dist = min(avgs_list)
        #print(best_dist)
        for num, dist in enumerate(avgs_list):
            if dist == best_dist:
                #print(taxa_list)
                #print(avgs_list)
                best_ref = taxa_list[num]
                return best_ref

def ref_selector(splits, phylo_distance_matrix, taxa_list):
    best_refs_list = []
    sorted_taxa = taxa_list[::-1]
    num_taxa = len(sorted_taxa)
    for split in splits:
        #print("NEWSPLIT")
        included_taxa = []
        num_ones = split.count('1')
        num_zeroes = split.count('0')
        if num_ones <= num_zeroes:
            for loc, position in enumerate(split):
                if position == '1':
                    #print(sorted_taxa[loc])
                    included_taxa.append(sorted_taxa[loc])
        elif num_ones > num_zeroes:
            for loc, position in enumerate(split):
                if position == '0':
                    #print(sorted_taxa[loc])
                    included_taxa.append(sorted_taxa[loc])

        best_ref = get_dists(included_taxa, phylo_distance_matrix)
        best_refs_list.append(best_ref)
    return best_refs_list


def main():
    args = parse_args()
    taxa = dendropy.TaxonNamespace()
    mle = dendropy.Tree.get(path=args.tree, schema='newick', taxon_namespace=taxa, preserve_underscores=True)
    mle_len = mle.length()
    mle.encode_bipartitions()
    pdc = mle.phylogenetic_distance_matrix()
    
    #for i, t1 in enumerate(mle.taxon_namespace[:-1]):
    #    for t2 in mle.taxon_namespace[i+1:]:
    #        print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, pdc(t1, t2)))

    tax_list = [] 
    #print(taxa)
    for taxon in taxa:
        tax_list.append(taxon)

    second_splits_branch_length = 0.0
    second_splits_branch_length_percent_of_tree = 0.0
    split_len_ratios = []
    grouped_splits = []
    split_n_lens = {}
    total_taxa = 0
    branch_mean = 0.0
    branch_std = 0.0
    nine_five_cutoff = 0.0
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
        split_len_ratios.append(split_edge_ratio)
        split_n_lens[str_bipart] = split_edge_ratio
        grouped_splits.append(str_bipart)
        for taxa in str_bipart:
            #print(taxa)
            taxa_in_split_count+=1
        total_taxa = taxa_in_split_count
    if args.dist == 'short': 
        branch_mean = numpy.mean(split_len_ratios)
    
        branch_std = numpy.std(split_len_ratios)
    
        nine_five_cutoff = branch_mean + (branch_std)
    elif args.dist == 'long':
        branch_mean = numpy.mean(split_len_ratios)

        branch_std = numpy.std(split_len_ratios)

        nine_five_cutoff = branch_mean + (branch_std + branch_std)

    #TAKES SPLITS AND ORGANIZES THEM BY HOW MANY TAXA THEY CONTAIN
    sort_splits = split_organiser(split_n_lens)
    #print(sort_splits)

    index = split_indexer(sort_splits)
    #print(index)

    nest_splits = nested_split_constructor(index, sort_splits, split_n_lens, nine_five_cutoff)
    #print(nest_splits)
    #print(tax_list)

    pick_refs = ref_selector(nest_splits, pdc, tax_list)
    for ref in pick_refs:
        print(ref)

    


if __name__ == '__main__':
    main()
