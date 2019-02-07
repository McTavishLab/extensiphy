#! /usr/bin/python
# script to import the output trees from the gon_phyling and phycorder trimmed
# and compare those two trees
# this script will also import pre-computed trees
# to assert you are computing roughly the same trees each time the test is run
# currently based on example tree-compare script from dendropy
# all thanks to the dendropy devs

import os
import argparse
import dendropy
from dendropy.calculate import treecompare

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--working_dir')
    parser.add_argument('--gon_phy_tree')
    parser.add_argument('--phycorder_tree')
    parser.add_argument('--orig_gon_phy_tree')
    parser.add_argument('--orig_phycorder_tree')
    return parser.parse_args()

def main():
    args = parse_args()

    # change into working dir where all the files are
    os.chdir(args.working_dir)

    # import the new gon_phyling produced tree
    gon_phy_tree = dendropy.Tree.get(
        path=args.gon_phy_tree,
        schema='newick',
        terminating_semicolon_required=False)

    # import the new phycorder produced tree
    phycorder_tree = dendropy.Tree.get(
        path=args.phycorder_tree,
        schema='newick',
        terminating_semicolon_required=False)

    # import the original gon_phyling produced tree
    orig_gon_phy_tree = dendropy.Tree.get(
        path=args.orig_gon_phy_tree,
        schema='newick',
        terminating_semicolon_required=False)

    # import the new phycorder produced tree
    orig_phycorder_tree = dendropy.Tree.get(
        path=args.orig_phycorder_tree,
        schema='newick',
        terminating_semicolon_required=False)

    #print(phycorder_tree)

    # establish common taxon namespace
    tns = dendropy.TaxonNamespace()

    ## ensure all trees loaded use common namespace
    tree1 = dendropy.Tree.get(
            data=gon_phy_tree,
            schema='newick',
            taxon_namespace=tns,
            terminating_semicolon_required=False)

    tree2 = dendropy.Tree.get(
            data=phycorder_tree,
            schema='newick',
            taxon_namespace=tns,
            terminating_semicolon_required=False)

    orig_tree1 = dendropy.Tree.get(
            data=orig_gon_phy_tree,
            schema='newick',
            taxon_namespace=tns,
            terminating_semicolon_required=False)

    orig_tree2 = dendropy.Tree.get(
            data=orig_phycorder_tree,
            schema='newick',
            taxon_namespace=tns,
            terminating_semicolon_required=False)

    # unweighted robinson-Foulds distance of the original gon_phyling (traditional)
    # tree compared to the newly produced gon_phyling tree
    print('\n')
    print("UNWEIGHTED RF distance comparison between majority rule consensus gon_phyling tree and the newly produced tree")
    print("RF: ")
    print(treecompare.symmetric_difference(orig_tree1, tree1))

    # unweighted robinson-Foulds distance of the original phycorder (rapid-updating)
    # tree compared to the newly produced phycorder tree
    print('\n')
    print("UNWEIGHTED RF distance comparison between majority rule consensus phycorder tree and the newly produced tree")
    print("RF: ")
    print(treecompare.symmetric_difference(orig_tree2, tree2))

    # Unweighted Robinson-Foulds distance
    print('\n')
    print("UNWEIGHTED RF distance comparison between majority rule consensus rapid-updating method (phycorder)")
    print("and traditional phylogenetics method: ")
    print("RF: ")
    majority_rule = treecompare.symmetric_difference(tree1, tree2)
    print(majority_rule)
    assert (majority_rule == int(2))



if __name__ == '__main__':
    main()
