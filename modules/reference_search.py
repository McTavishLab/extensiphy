#! /usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import statistics
import numpy
import dendropy

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', help='phylogenetic tree file in newick format')
    parser.add_argument('--dist', default='short', help='branch length distance leading to references (DEFAULT:short)(Options: short, long')
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

#MAIN FUNCTION TO ITERATE THROUGH SPLITS AND DECIDE WHAT SPLITS ARE INCLUDED FOR FINDING REFERENCES
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
<<<<<<< HEAD
    
    #Test tree for dev purposes.
    mle_1 = '(taxon_10:0.000490504536564,(taxon_15:0.00147793977114,((taxon_22:0.00775296610854,(taxon_11:0.00233219469227,taxon_19:0.00292833772938):0.00661072260997):0.00437388771198,((taxon_17:0.00238143585578,((taxon_20:0.00040300210567,((taxon_27:0.000706168118649,taxon_23:0.000100735490224):1.00000050003e-06,(taxon_25:0.000503957096912,(taxon_26:1.00000050003e-06,taxon_21:0.000100729540367):0.000403230214661):1.00000050003e-06):1.00000050003e-06):0.000200876170437,(taxon_28:0.0011141033398,taxon_24:9.95483145232e-05):0.000507160142818):0.00184911046943):0.00339928596404,(taxon_16.ref:0.000400637478349,(taxon_13:0.000100794665489,(taxon_14:0.000706692317158,(taxon_18:1.00000050003e-06,taxon_1:1.00000050003e-06):0.000302521801106):1.00000050003e-06):0.00061029754872):0.00511761213428):0.000906053315062):0.0040106698136):0.00276938494406,taxon_12:0.00214526732704):0.0;'
    
    #Establish namespace.
    taxa = dendropy.TaxonNamespace()
    
    mle = ''

    #Handle whether working with test tree or actual input tree.
    if args.test == "True":
        mle = dendropy.Tree.get(data=mle_1, schema='newick', taxon_namespace=taxa)
    else:
        mle = dendropy.Tree.get(path=args.tree_file, schema='newick', taxon_namespace=taxa)
    

    mle_len = mle.length()
    print("tree length is", mle_len)
    for split in mle.encode_bipartitions():
        print(split.leafset_bitmask)
        print(dir(split))
        for split2 in mle.encode_bipartitions():
            if split.is_nested_within(split2):
                print(split.leafset_as_newick_string(taxa))
                print(split2.leafset_as_newick_string(taxa))
                print("is nested")
    

    #print(dir(mle))
    
    
    # taxa_count = 0
    # for name in taxa:
    #     taxa_count+=1

    # split_count = 0
    # for split in mle.bipartition_encoding:
    #     split_count+=1

    # #taxa = dendropy.TaxonNamespace()
    # new_mle = dendropy.Tree.get(data=mle_1, schema='newick', taxon_namespace=taxa)
    # new_mle.encode_bipartitions()
    # split_count_2 = 0
    # for split in mle.bipartition_encoding:
    #     split_count_2+=1
    #     if split_count_2 == split_count:
    #         print(dir(mle.bipartition_edge_map[split]))
    #     #print(split.edge)
    #     #print(split)
    #     print(mle.bipartition_edge_map[split].length)
    #     #print(dir(mle.bipartition_edge_map[split]))
    #     print(mle.bipartition_edge_map[split].leafset_bitmask)
    #     print(mle.bipartition_edge_map[split].bipartition)
    #     #print(split.split_bitmask)
    #     #print(split.split_as_newick_string.__getattribute__)
    #     #print(dir(split.split_as_newick_string))
    #     #print(split.edge)
    #     #print(type(split))
    #     #print(dir(split))
    #     pass

    # #print(mle.bipartition_edge_map)
    # #for edge in mle.bipartition_edge_map:
    # #    print(edge)

    # print(taxa)
    # #print(type(taxa))
    # #print(dir(taxa))
    # #for name in taxa:
    # #    print(dir(name))

if __name__ == '__main__':
    main()
