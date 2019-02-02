#! /usr/bin/python
# script to import the output trees from the gon_phyling and phycorder trimmed
# and compare those two trees
# this script will also import pre-computed trees
# to assert you are computing roughly the same trees each time the test is run

import dendropy
from dendropy.calculate import treecompare

s1 = "(a,(b,(c,d)));"
s2 = "(a,(d,(b,c)));"

# establish common taxon namespace
tns = dendropy.TaxonNamespace()

# ensure all trees loaded use common namespace
tree1 = dendropy.Tree.get(
        data=s1,
        schema='newick',
        taxon_namespace=tns)
tree2 = dendropy.Tree.get(
        data=s2,
        schema='newick',
        taxon_namespace=tns)

## Unweighted Robinson-Foulds distance
print(treecompare.symmetric_difference(tree1, tree2))
