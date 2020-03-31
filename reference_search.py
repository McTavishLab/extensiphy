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

#SEPARATE OUT ALL SPLITS CONTAINING A SINGLE TAXON AND RETURN THAT LIST OF SPLITS
def separate_single_taxon_splits(all_biparts):
    single_taxon_node_list = []
    for split in all_biparts:
        taxon_count = 0
        for taxon in split:
            if taxon == '1':
                taxon_count+=1
        if taxon_count == 1:
            single_taxon_node_list.append(split)
    return single_taxon_node_list

def separate_remaining_splits(all_biparts, large_splits, single_taxon_splits):
    remaining_splits = []
    for split in all_biparts:
        if split not in large_splits:
            if split not in single_taxon_splits:
                remaining_splits.append(split)
    return remaining_splits

def combine_splits(split_1, split_2):
    list_of_combined_taxa_positions = list(map(list, zip(split_1, split_2)))
    return list_of_combined_taxa_positions


def split_nesting(filtered_split_list):
    big_split_list = []
    current_split = ''
    previous_split = ''
    first_split = ''
    split_count = 0
    first_split = 0
    nested_count = 0
    previous_split_count = 0
    for split in filtered_split_list:
                
        if first_split == 0:
            first_split+=1
            first_split = split
            #previous_split = split
        elif first_split >= 1:
            #print(filtered_split_list[split_count])
            #print(split)
            prev_split_loc = int(split_count - 1)
            #print(prev_split_loc)
            #print(split_count)
            current_split = split
            previous_split = filtered_split_list[split_count]
            #combine = list(map(list, zip(current_split, previous_split)))
            combine = combine_splits(current_split, previous_split)
            #print(combine)
            check_combine = list(map(same_taxa, combine))
            #print(''.join(check_combine))
            #print(check_combine)
            #print(list(map(type, check_combine)))
            int_prev = list(map(int, previous_split))
            #print(int_prev)
            if check_combine == int_prev:
                #DO SOMETHING WITH SPLIT BECAUSE ITS NESTED
                #print("nested")
                nested_count+=1
            else:
                #print("not nested")
                current_nested_splits_list = []
                for specific_split in range((split_count - nested_count), (split_count + 1)):
                    #print(filtered_split_list[specific_split])
                    current_nested_splits_list.append(filtered_split_list[specific_split])
                    #TODO: HANDLE MULTIPLE INNER NESTED SETS OF SPLITS NOT STACKED ON TOP OF EACHOTHER
                nested_count = 0
                big_split_list.append(current_nested_splits_list)
            #print("reset")
            split_count+=1
    return big_split_list

def same_taxa(taxon_positions_list):
    if taxon_positions_list[1] == '1':
        if taxon_positions_list[0] == taxon_positions_list[1]:
            #print("SAME")
            return 1
        else:
            return 0
    else:
        #print("DIFF")
        return 0

#TAKE A SET OF NESTED SPLITS THAT DONT CONTAIN THE SPLITS LEADING TO A SINGLE TAXON AND A LIST OF SINGLE TAXON SPLITS
#MATCH SINGLE TAXON SPLITS TO APPROPRIATE NESTED SPLITS
def nest_leaves(single_taxon_splits, nested_splits):
    for nested_set in nested_splits:
        for bigger_split in nested_set:
            for split in single_taxon_splits:
                combine = list(map(list, zip(bigger_split, split)))
                check_combine = list(map(same_taxa, combine))
                #print(check_combine)
                if 1 in check_combine:
                    if split not in nested_set:
                        nested_set.append(split)
                else:
                    continue
    return(nested_splits)
    #for nested_set in nested_splits:
    #    #print("##########################################################################")
    #    for split in nested_set:
    #        print(split)

#TODO: KEEP TRACK OF NESTED SPLITS IN RELATION TO EACHOTHER
# TO KNOW WHICH SPLITS FIT INTO EACHOTHER
def nested_lens(splits_with_lens_dict, nested_splits_list, tree_len):
    nested_splits_with_len_ratio_list = []
    #print(splits_with_lens_dict)
    #print(nested_splits_list)
    print(tree_len)
    for num, split_set in enumerate(nested_splits_list):
        splits_with_len_ratio_dict = {}
        #print(split_set)
        for num_2, split in enumerate(split_set):
            #split_n_len_ratio = []
            #print(splits_with_lens_dict[split])
            split_edge_ratio = float(splits_with_lens_dict[split]) / float(tree_len)
            #split_n_len_ratio.append(split)
            #splitn_len_ratio.append(split)
            splits_with_len_ratio_dict[split] = split_edge_ratio
        #print("end of nested split")
        nested_splits_with_len_ratio_list.append(splits_with_len_ratio_dict)
    return nested_splits_with_len_ratio_list

def split_type_sorter(dict_of_splits_n_len_ratios, group_count):
    output_dict = {}
    multiple_taxa_long_branch_split_count = 0
    multiple_taxa_long_branch_split_list = []
    single_taxa_long_branch_count = 0
    single_taxa_long_branch_list = []
    no_long_branches_count = 0
    no_long_branches_list = []
    for split, len_ratio in dict_of_splits_n_len_ratios.items():
        taxa_count = 0
        for position in split:
            if position == '1':
                taxa_count+=1
        if taxa_count > 1:
            if float(len_ratio) >= float(1.0):
                multiple_taxa_long_branch_split_count+=1
                multiple_taxa_long_branch_split_list.append(split)
            elif float(len_ratio) > float(1.0):
                no_long_branches_count+=1
                no_long_branches_list.append(split)
        elif taxa_count == 1:
            if float(len_ratio) >= float(1.0):
                single_taxa_long_branch_count+=1
                single_taxa_long_branch_list.append(split)
            elif float(len_ratio) < float(1.0):
                no_long_branches_count+=1
                no_long_branches_list.append(split)

    if multiple_taxa_long_branch_split_count == 1:
        if single_taxa_long_branch_count == 0:
            output_dict['status'] = 'long_branch_to_multiple_short'
            output_dict['short_branches'] = no_long_branches_list
            output_dict['long_branch'] = multiple_taxa_long_brach_split_list




#TAKE NESTED SETS OF SPLITS WITH BRANCH LENGTH RATIOS
#OUTPUT BEST SPLITS BASED ON LONG BRANCHES AND SPLIT SETS THAT NEED REFERENCES CHOSEN BY DISTANCES
#THIS IS UGLY A HELL
def output_refs(list_of_nested_splits_len_ratios):
    output_dict = {}
    for group in list_of_nested_splits_len_ratios:
        print("group")
        groups_split_info_dict = {"sngle_tax_lng_brnch":[], "shrt_brnch_tax":[]}
        multiple_taxa_long_branch_split_count = 0
        for split, len_ratio in group.items():
            tax_in_split_count = 0
            print("split '%s' and and branch length ratio '%s'" % (split, len_ratio))
            for position in split:
                if position == '1':
                    tax_in_split_count+=1
            if tax_in_split_count > 1:



        


#SEPARATE OUT ALL SPLITS CONTAINING 55% OR MORE OF TAXA,
#THESE SPLITS CONTAIN SO MANY TAXA THAT BRANCH LENGTHS MUST BE SHORT 
#FOR A REFERENCE TO REPRESENT ALL TAXA IN THE SPLIT
def separate_large_splits(all_biparts, total_taxa_num):
    #print(total_taxa_num)
    large_split_list = []
    cutoff = 0.55
    for split in all_biparts:
        taxon_count = 0
        #print(type(split[0]))
        #print(taxon_count)
        in_taxon = "1"
        for taxon in split:
            #print(in_taxon)
            #print(type(in_taxon))
            #print(taxon)
            #print(type(taxon))
            #print(int(taxon) + int(in_taxon))
            if int(taxon) == int(in_taxon):
                taxon_count+=1
                #print("WAFFLE")
        #print(taxon_count)
        proportion_of_all_taxa = float(float(taxon_count) / float(total_taxa_num))
        #print(proportion_of_all_taxa)
        if proportion_of_all_taxa > cutoff:
            large_split_list.append(split)
    return large_split_list

def count_taxa(namespace_list):
    total_taxa_count = 0
    #print(namespace_list)
    for taxon in namespace_list:
        total_taxa_count+=1
        print(total_taxa_count)
    return total_taxa_count

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
    for i, t1 in enumerate(mle.taxon_namespace[:-1]):
        for t2 in mle.taxon_namespace[i+1:]:
            print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, pdc(t1, t2)))

    
    print(taxa)


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
    #for split in multiple_taxon_leaf_splits:
    #    print(split)
    #pair_bipartitions = splits_iter(split_n_lens, grouped_splits)
   

    #print(total_taxa)
    tax_count = count_taxa(taxa) 
    #print(tax_count)
    
    #print("WAFFLE")
    big_splits = separate_large_splits(grouped_splits, total_taxa)
    #print(big_splits)

    single_taxons = separate_single_taxon_splits(grouped_splits)

    workable_splits = separate_remaining_splits(grouped_splits, big_splits, single_taxons) 
    #for split in workable_splits:
    #    print(split)
    #print("WAFFLE2")
    
    split_analysis = split_nesting(workable_splits)
    #for split in split_analysis:
    #    print(split)

    #ADD BACK IN SINGLE TAXON SPLITS AND FIND WHICH NEST GROUP THEY BELONG TO
    add_singles = nest_leaves(single_taxons, split_analysis)

    find_nested_lens = nested_lens(split_n_lens, add_singles, mle_len)
    for nested_splits in find_nested_lens:
        print(nested_splits)

    get_refs = output_refs(find_nested_lens)





if __name__ == '__main__':
    main()
