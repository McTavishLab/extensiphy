import dendropy
import sys

mat=sys.argv[1]
mattype=sys.argv[2]
tre=sys.argv[3]
tretype=sys.argv[4]
nam=sys.arg[5]


charmat = dendropy.DataSet.get_from_path(fas, mattype)
charmat.write_to_path("{}.fas".format(nam), "dnafasta")
charmat.write_to_path("{}.phy".format(nam), "phylip")


tre=dendropy.Tree.get_from_path(tre,tretype)
tre.resolve_polytomies()
tre.write_to_path("{}.tre".format(nam),"newick",quote_underscores=False)
