#! /usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import statistics
import numpy
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
    #sorted_index = index[::-1]
    index.sort()
    #return sorted_index
    return index


def nested_split_constructor(index, organized_splits, splits_and_ratios_dict, cutoff):
    final_splits_list = []
    #print(cutoff)
    for num in index:
        if num != 0:
            #print(num)
            for split in organized_splits[num]:
                #print(splits_and_ratios_dict[split])
                lookat_split = splits_and_ratios_dict[split]
                if lookat_split >= cutoff:
                    #print("found")
                    final_splits_list.append(split)
                    
    return final_splits_list


def get_dists(taxa_list, phylo_distance_matrix):
    finish = 0
    dist_list = []
    avgs_list = []
    single_taxa = []
    best_dist = 0
    best_ref = ''
    for subject_tax in taxa_list:
        for selected_tax in taxa_list:
            if selected_tax != subject_tax:
                dists = phylo_distance_matrix(subject_tax, selected_tax)
                dist_list.append(dists)
                avg = numpy.mean(dist_list)
                avgs_list.append(avg)
            elif len(taxa_list) == 1:
                best_ref = subject_tax
                finish = 1
                return best_ref
    if finish != 1:
        best_dist = min(avgs_list)
        for num, dist in enumerate(avgs_list):
            if dist == best_dist:
                best_ref = taxa_list[num]
                return best_ref

def ref_selector(splits, phylo_distance_matrix, taxa_list):
    sorted_taxa = taxa_list[::-1]
    num_taxa = len(sorted_taxa)
    for split in splits:
        print("NEWSPLIT")
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

        best_dists = get_dists(included_taxa, phylo_distance_matrix)
        print(best_dists)


def main():
    args = parse_args()
    taxa = dendropy.TaxonNamespace()
    mle = dendropy.Tree.get(path=args.tree_file, schema='newick', taxon_namespace=taxa, preserve_underscores=True)
    mle_len = mle.length()
    mle.encode_bipartitions()
    pdc = mle.phylogenetic_distance_matrix()
    
    #for i, t1 in enumerate(mle.taxon_namespace[:-1]):
    #    for t2 in mle.taxon_namespace[i+1:]:
    #        print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, pdc(t1, t2)))

    tax_list = [] 
    print(taxa)
    for taxon in taxa:
        tax_list.append(taxon)

    second_splits_branch_length = 0.0
    second_splits_branch_length_percent_of_tree = 0.0
    split_len_ratios = []
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
        split_len_ratios.append(split_edge_ratio)
        split_n_lens[str_bipart] = split_edge_ratio
        grouped_splits.append(str_bipart)
        for taxa in str_bipart:
            #print(taxa)
            taxa_in_split_count+=1
        total_taxa = taxa_in_split_count
    #print(split_len_ratios)   
    branch_mean = numpy.mean(split_len_ratios)
    #print(branch_mean)
    branch_std = numpy.std(split_len_ratios)
    #print(branch_std)
    nine_five_cutoff = branch_mean + (branch_std)
    print(nine_five_cutoff)
    sort_splits = split_organiser(split_n_lens)
    #print(sort_splits)

    index = split_indexer(sort_splits)
    #print(index)

    nest_splits = nested_split_constructor(index, sort_splits, split_n_lens, nine_five_cutoff)
    #print(nest_splits)
    #print(tax_list)

    pick_refs = ref_selector(nest_splits, pdc, tax_list)
    print(pick_refs)

    


if __name__ == '__main__':
    main()
