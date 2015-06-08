import dendropy
import sys
"""resolves polytomies in tree, writes character dataset out to fasta"""

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
        print("Taxon '{}': found in both trees".format(taxon.label))