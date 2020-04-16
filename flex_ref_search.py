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

def position_check(list_of_2_taxa):
    if list_of_2_taxa[0] == '1':
        if list_of_2_taxa[0] == list_of_2_taxa[1]:
            return 1
        else:
            return 0
    else:
        return 0

def combine_splits(split_1, split_2):
    list_of_combined_taxa_positions = list(map(list, zip(split_1, split_2)))
    which_nests = list(map(position_check, list_of_combined_taxa_positions))
    which_nests = ''.join(str(e) for e in which_nests)

    return which_nests

def single_tax_long_branch_finder(index, organized_splits, splits_and_ratios_dict, cutoff):
    final_splits_list = []
    #print(cutoff)
    for num in index:
        if num == 1:
            #print(num)
            for split in organized_splits[num]:
                #print(splits_and_ratios_dict[split])
                lookat_split = splits_and_ratios_dict[split]
                if lookat_split >= cutoff:
                    #print("found")
                    final_splits_list.append(split)

    return final_splits_list

#MAIN FUNCTION TO ITERATE THROUGH SPLITS AND DECIDE WHAT SPLITS ARE INCLUDED FOR FINDING REFERENCES
def nested_split_constructor(index, organized_splits, splits_and_ratios_dict, cutoff):
    final_splits_list = []
    paths = {}
    for split in organized_splits[1]:
        split_path = {}
        #print("START HERE")
        #print(split)
        current_split = split
        next_split = ''
        for num in index[2:]:
            #print(num)
            for bigger_split in organized_splits[num]:
                next_split = bigger_split
                nest_check = combine_splits(current_split, next_split)
                if nest_check == current_split:
                    #print(current_split)
                    #print(next_split)
                    split_path[num] = next_split
        paths[split] = split_path
    
    splits_to_examine = []
    for split, path in paths.items():
        splits_passed_filter = []
        #print("WAFFLE")
        #print(split)
        #for key, value in path.items():
        #    print(key)
        #    print(value)
        aggregate_branch = 0.0
        splits_switch = 0
        for split_level, nested_split in path.items():
            num_taxa = 0
            #splits_switch = 0
            for position in nested_split:
                if position == '1':
                    num_taxa+=1
            current_branch_len = splits_and_ratios_dict[nested_split]
            aggregate_branch = aggregate_branch + current_branch_len
            #print(aggregate_branch)
            if aggregate_branch >= cutoff:
                if splits_switch == 0:
                    #print("LEN CUTOFF")
                    splits_switch = 1
                    #print(nested_split)
                    if nested_split not in splits_to_examine:
                        splits_to_examine.append(nested_split)
            elif (num_taxa / max(index)) >= 0.15:
                if splits_switch == 0:
                    #print("num taxa cutoff")
                    splits_switch = 1
                    #print(nested_split)
                    if nested_split not in splits_to_examine:
                        splits_to_examine.append(nested_split)
    #print(splits_to_examine)



    return splits_to_examine


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




#def spread_refs(spread_num, ref_list, organized_split_dict, phylo_distance_matrix, taxa_list, splits_producing_refs):
#    sorted_taxa = taxa_list[::-1]
#    all_split_levels = []
#    for key, value in organized_split_dict.items():
#        all_split_levels.append(key)
#    #print(all_split_levels)
#    closest_split = min(all_split_levels, key=lambda x:abs(x-spread_num))
#    
#    splits_to_keep = {}
#    analyzed_taxa = {}
#    positions_used = []
#    for split in splits_producing_refs:
#        taxa_from_current_split = []
#        count_tax = split.count('1')
#        if count_tax <= closest_split:
#            for pos, taxon in enumerate(split):
#                if taxon == '1':
#                    taxa_from_current_split.append(taxa_list[pos])
#                    if pos not in positions_used:
#                        positions_used.append(pos)
#        analyzed_taxa[split] = taxa_from_current_split
#    #for split, tax_set in analyzed_taxa.items():
#    #    print(split)
#    #    print(tax_set)
#    print(positions_used)
#    tax_groups = []
#    group = []
#    for num, tax_position in enumerate(positions_used):
#        if num != 0:
#            if tax_position == positions_used[num - 1] + 1:
#                print(tax_position)
#            else:
#                print("new group")
        


def spread_refs(refs, phylo_dist_matrix):
    dist_dict = {}
    total_dists = []
    for i, t1 in enumerate(refs[:-1]):
        dists_from_t1 = []
        for t2 in refs[i+1:]:
            dist_from_dict = {}
            #print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, phylo_dist_matrix(t1, t2)))
            dist_from_dict[t2] = phylo_dist_matrix(t1, t2)
            total_dists.append(phylo_dist_matrix(t1, t2))
            dists_from_t1.append(dist_from_dict)
        dist_dict[t1] = dists_from_t1

    ref_dist_mean = numpy.mean(total_dists)
    ref_dist_std = numpy.std(total_dists)

    #print(ref_dist_mean)
    #print(ref_dist_std)
    tax_to_keep = []
    for main_tax, tax_set in dist_dict.items():
        for tax_dist_pair in tax_set:
            for key, value in tax_dist_pair.items():
                if value <= (ref_dist_mean - ref_dist_std):
                    #print("to close")
                    pass
                elif value > (ref_dist_mean - ref_dist_std) and key not in tax_to_keep:
                    tax_to_keep.append(key)
    #print(tax_to_keep)

    return tax_to_keep

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
    #for num, split_set in sort_splits.items():
    #    print(num)
    #    print(split_set)

    index = split_indexer(sort_splits)
    #print(index)

    nest_splits = nested_split_constructor(index, sort_splits, split_n_lens, nine_five_cutoff)
    #print(nest_splits)
    #print(tax_list)

    single_tax_branches = single_tax_long_branch_finder(index, sort_splits, split_n_lens, nine_five_cutoff)    

    nest_pick_refs = ref_selector(nest_splits, pdc, tax_list)
    nest_pick_refs = set(nest_pick_refs)
    nest_pick_refs_list = list(nest_pick_refs)
    #for ref in nest_pick_refs:
    #    print(ref)
    #print(nest_pick_refs)
    #print("SINGLE TAX REFS HERE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")

    single_tax_refs = ref_selector(single_tax_branches, pdc, tax_list)
    single_tax_refs = set(single_tax_refs)
    #for ref in single_tax_refs:
    #    print(ref)

    unique_refs = nest_pick_refs.union(single_tax_refs)
    list_refs = list(unique_refs)
    print(list_refs)
    for ref in list_refs:
        print(ref)

    #for i, t1 in enumerate(list_refs[:-1]):
    #    for t2 in list_refs[i+1:]:
    #        print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, pdc(t1, t2)))

    #gold_spread = int((max(index) * 0.15))
    
    #spread_out_refs = spread_refs(gold_spread, list_refs, sort_splits, pdc, tax_list, nest_splits)

    trim_refs = spread_refs(nest_pick_refs_list, pdc)
    #print(trim_refs)
    #for ref in trim_refs:
    #    print(ref)
if __name__ == '__main__':
    main()
