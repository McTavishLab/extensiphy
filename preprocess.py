#! /usr/bin/env python
from dendropy import DnaCharacterMatrix, Tree
import sys

mat=sys.argv[1]
mattype=sys.argv[2]
tre=sys.argv[3]
tretype=sys.argv[4]
nam=sys.arg[5]

mat = 'example.aln'
mattype = 'fasta'
tre = 'tree.tre'
tretype = 'newick'

d = DnaCharacterMatrix.get(path=mat,
                           schema=mattype)
# make the taxon_namespace immutable, so the tree does not add
#   new labels...
d.taxon_namespace.is_mutable = False
tree = Tree.get(path=tre,
                schema=tretype,
                preserve_underscores=True,
                taxon_namespace=d.taxon_namespace)

# get all of the taxa associated with tips of the tree, and make sure that
#   they include all of the members of the data's taxon_namespace...
treed_taxa = [i.taxon for i in tree.leaf_nodes()]
if len(treed_taxa) != len(d.taxon_namespace):
    missing = [i.label for i in d.taxon_namespace if i not in treed_taxa]
    emf = 'Some of the taxa are not in the tree. Missing "{}"\n'
    em = emf.format('", "'.join(missing))
    raise ValueError(em)


"""resolves polytomies in tree, writes character dataset out to fasta

mat=sys.argv[1]
mattype=sys.argv[2]
tre=sys.argv[3]
tretype=sys.argv[4]
nam=sys.arg[5]

mat="example.aln"
mattype="fasta"
tre="tree.tre"
tretype="newick"
nam="test"
d = dendropy.DnaCharacterMatrix.get(path=mat, schema=mattype)
tree = dendropy.Tree.get(path='tree.tre', schema='newick', preserve_underscores=True, taxon_namespace=d.taxon_namespace)


tre.resolve_polytomies()
tre.write_to_path("{}.tre".format(nam), "newick", quote_underscores=False, )


print(d.taxon_namespaces[0].description(2))


charmat.write_to_path("{}.fas".format(nam), ")
#charmat.write_to_path("{}.phy".format(nam), "phylip")






for taxon in charmat.taxon_namespace:
    if taxon in tree_list2.taxon_namespace:
        # this branch is never visited
        print("Taxon '{}': found in both trees".format(taxon.label))"""
