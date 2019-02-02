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
    return parser.parse_args()

def main():
    args = parse_args()

    # change into working dir where all the files are
    os.chdir(args.working_dir)

    open_gon_phy = open(args.gon_phy_tree, 'r')
    read_gon_phy = open_gon_phy.read()
    #print(read_gon_phy)

    open_phycord = open(args.phycorder_tree, 'r')
    read_phycord = open_phycord.read()
    #print(read_phycord)

    # import gon_phyling produced tree
    gon_phy_tree = dendropy.Tree.get(
        path=args.gon_phy_tree,
        schema='newick',
        terminating_semicolon_required=False)

    #print(gon_phy_tree)

    # import phycorder produced tree
    phycorder_tree = dendropy.Tree.get(
        path=args.phycorder_tree,
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
    #
    # Unweighted Robinson-Foulds distance
    print('\n')
    print("UNWEIGHTED RF distance: ")
    print(treecompare.symmetric_difference(tree1, tree2))


if __name__ == '__main__':
    main()
