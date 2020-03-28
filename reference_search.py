#! /usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import dendropy

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree_file')
    return parser.parse_args()


def splits_iter(splits_and_len_dict, bipart_list):
    #LOOP THROUGH DICT OF SPLITS AND BRANCH LENGTHS AS PERCENTAGES OF TOTAL TREE LENGTH
    #CREATE LIST OF TUPLES CONTAINING TWO BIPARTITIONS AS LONG AS THE BIPARTITIONS ARE NOT IDENTICAL
    tuple_list = []
    for bipart in bipart_list:
        #print(bipart)
        for split, len_ratio in splits_and_len_dict.items():
            if bipart != split:
                bipart_tuple = (bipart, split)
                tuple_list.append(bipart_tuple)
    return tuple_list

def splits_compare(split_tuple_list, split_and_branch_len_dict):
    #COMPARE BIPARTITIONS TO IDENTIFY NESTED BIPARTITIONS
    for bipart_tuple in split_tuple_list:
        first_bipart_positions = []
        second_bipart_positions = []

        #for pos, value in enumerate(bipart_tuple[0]):
        #    first_bipart_positions[pos] = value
        #for pos, value in enumerate(bipart_tuple[1]):
        #    second_bipart_positions[pos] = value
        #print(first_bipart_positions)
        #print(second_bipart_positions)
        #print("WAFFLE")
        

def main():
    args = parse_args()
    mle_1 = '(taxon_10:0.000490504536564,(taxon_15:0.00147793977114,((taxon_22:0.00775296610854,(taxon_11:0.00233219469227,taxon_19:0.00292833772938):0.00661072260997):0.00437388771198,((taxon_17:0.00238143585578,((taxon_20:0.00040300210567,((taxon_27:0.000706168118649,taxon_23:0.000100735490224):1.00000050003e-06,(taxon_25:0.000503957096912,(taxon_26:1.00000050003e-06,taxon_21:0.000100729540367):0.000403230214661):1.00000050003e-06):1.00000050003e-06):0.000200876170437,(taxon_28:0.0011141033398,taxon_24:9.95483145232e-05):0.000507160142818):0.00184911046943):0.00339928596404,(taxon_16.ref:0.000400637478349,(taxon_13:0.000100794665489,(taxon_14:0.000706692317158,(taxon_18:1.00000050003e-06,taxon_1:1.00000050003e-06):0.000302521801106):1.00000050003e-06):0.00061029754872):0.00511761213428):0.000906053315062):0.0040106698136):0.00276938494406,taxon_12:0.00214526732704):0.0;'
    taxa = dendropy.TaxonNamespace()
    mle = dendropy.Tree.get(path=args.tree_file, schema='newick', taxon_namespace=taxa)
    mle_len = mle.length()
    #print(dir(mle))
    #print(mle.bipartition_edge_map)
    #print("tree length is", mle_len)
        #print(dir(mle))
        #print(mle.bipartition_edge_map)
    mle.encode_bipartitions()
    pdc = mle.phylogenetic_distance_matrix()

    
    taxa_count = 0
    for name in taxa:
        taxa_count+=1

    split_count = 0
    for split in mle.bipartition_encoding:
        split_count+=1

    #taxa = dendropy.TaxonNamespace()
    new_mle = dendropy.Tree.get(data=mle_1, schema='newick', taxon_namespace=taxa)
    new_mle.encode_bipartitions()
    split_count_2 = 0
    splits_branch_length = 0.0
    splits_branch_length_percent_of_tree = 0.0
    for split in mle.bipartition_encoding:
        split_count_2+=1
        #if split_count_2 == split_count:
            #print(dir(mle.bipartition_edge_map[split]))
        #print(split.edge)
        #print(split)
        #print(mle.bipartition_edge_map[split].length)
            #print(dir(mle.bipartition_edge_map[split]))
            #print(mle.bipartition_edge_map[split].leafset_bitmask)
        #print(mle.bipartition_edge_map[split].bipartition)
        
        split_edge_ratio = float(mle.bipartition_edge_map[split].length) / float(mle_len)
        #REMEMBER THAT BIPARTITION SPLITS ARE BACKWARDS COMPARED TO THE LIST PRODUCED BY NAMESPACE
        #print("branch length as ratio of tree: %f" % (split_edge_ratio))
        #CURRENT IDEA: LOOK AT SPLIT, SEE WHAT THE DISTANCE IS TO THE SPLIT THAT CONTAINS THIS SPLIT
        #IF LONG: USE SPLIT
        #IF SHORT: CONSIDER CONTAINER SPLIT
        splits_branch_length+=mle.bipartition_edge_map[split].length
        splits_branch_length_percent_of_tree+=split_edge_ratio
            #print(split.split_bitmask)
            #print(split.split_as_newick_string.__getattribute__)
            #print(dir(split.split_as_newick_string))
            #print(split.edge)
            #print(type(split))
            #print(dir(split))
        pass
    
    #print(splits_branch_length)
    #print(mle_len)
    #print(splits_branch_length_percent_of_tree)
    #for i, t1 in enumerate(mle.taxon_namespace[:-1]):
    #    for t2 in mle.taxon_namespace[i+1:]:
    #        print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, pdc(t1, t2)))

    
    print(taxa)

    second_splits_branch_length = 0.0
    second_splits_branch_length_percent_of_tree = 0.0
    grouped_splits = []
    split_n_lens = {}
    split_encode = mle.bipartition_encoding
    for split in split_encode:
        split_branch_len = mle.bipartition_edge_map[split].length
        bipart = mle.bipartition_edge_map[split].bipartition
        split_edge_ratio = float(split_branch_len) / float(mle_len)
        second_splits_branch_length+=mle.bipartition_edge_map[split].length
        second_splits_branch_length_percent_of_tree+=split_edge_ratio
        #print(type(bipart))
        str_bipart = str(bipart)
        print(str_bipart)
        #print(type(str_bipart))
        split_n_lens[str_bipart] = split_edge_ratio
        grouped_splits.append(str_bipart)
        
    pair_bipartitions = splits_iter(split_n_lens, grouped_splits)
    #print(pair_bipartitions)

    #bipart_analysis = splits_compare(pair_bipartitions, split_n_lens)



if __name__ == '__main__':
    main()
