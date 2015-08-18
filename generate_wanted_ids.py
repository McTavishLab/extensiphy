#!/usr/bin/env python
import sys
from dendropy import Tree
from peyotl import gen_otu_dict, iter_node
from peyotl.manip import iter_trees
from peyotl.phylesystem.phylesystem_umbrella import Phylesystem
from peyotl.sugar import tree_of_life, taxonomy
from peyotl.nexson_syntax import extract_tree,  PhyloSchema
import re


phy = Phylesystem()

out = codecs.getwriter('utf-8')(sys.stdout)
study_id = "pg_2739"
tree_id = "tree6601"
print study_id
n = phy.return_study(study_id)[0]

#There has gort to be a better way to get the ingroup otus...
m = extract_tree(n, tree_id, PhyloSchema('newick', output_nexml2json = '1.2.1', content="tree", tip_label="ot:ottId"), subtree_id="ingroup")

otu_dict = gen_otu_dict(n)
ottids = []
for oid, o in otu_dict.items():
    try:
        ottid = o[u'^ot:ottId']
        if ("{}:".format(ottid) in m) or ("{})".format(ottid) in m) or ("{},".format(ottid) in m):
            ottids.append(ottid)        
        else:
            print(o)
    except:
        pass


mrca_node = tree_of_life.mrca(ott_ids=ottids, wrap_response=True)

#curl -X POST http://api.opentreeoflife.org/v2/taxonomy/taxon -H "content-type:application/json" -d '{"ott_id":225495, "list_terminal_descendants": "True"}'

named_node = tree_of_life.mrca(ott_ids=[mrca_node.nearest_taxon.ott_id], wrap_response=True)
newick = named_node.subtree_newick # NOTE: Excludes daughter taxa that  are considered non monophylo.


ott_to_ncbi = {}
fi =open("ott_ncbi")
for lin in fi:
    lii= lin.split(",")
    ott_to_ncbi[lii[0]]=lii[1]


matches = re.finditer('ott(\d*)', newick) 
wanted_ncbi_ids = []
unfound = []
for item in matches:
    try:
        wanted_ncbi_ids.append(ott_to_ncbi[item.group(1)])
    except:
        unfound.append(item.group(1))


#NOW: work through the alignment, blast each seq (EACH?)
#and make a list of genbacnk id's that match.

#OR ancestral states?
#blastn -remote -query fulltest3/cns.fa -db nr -max_target_seqs 1 -outfmt "6 qseqid staxids"




#Now get all descendants of that node in the synth tree. Or the taxonomy? better closest tax node above.


'''



info = taxonomy.taxon(ott_id= mrca_node.nearest_taxon.ott_id,
                          include_lineage=False,
                          list_terminal_descendants=True,
                          wrap_response=False)
named_node = tree_of_life.mrca(ott_ids=[mrca_node.nearest_taxon.ott_id], wrap_response=True)
newick = named_node.subtree_newick

names = re.split(r'[;,\s]\s*', newick)
names = newick.replace('(',',').replace(')',',').split(",")













DOMAIN = 'https://tree.opentreeoflife.org'
URL_PATH_FMT ='opentree/argus/otol.draft.22@{i:d}'
URL_FMT = DOMAIN + '/' + URL_PATH_FMT

url = URL_FMT.format(i=mrca_node.node_id)        
if mrca_node.is_taxon:
        errstream.write('The node in the Graph of Life corresponds to a taxon:\n')
        mrca_node.write_report(errstream)
else:
    errstream.write('The node in the Graph of Life does not correspond to a taxon.\nThe most recent ancestor which is also a named taxon in OTT is:\n')
    mrca_node.nearest_taxon.write_report(errstream)
if subtree:
        # We could ask for this using: 
        #   newick = tree_of_life.subtree(node_id=mrca_node.node_id)['newick']
        # or we can ask the GoLNodeWrapper object to do the call (as shown below)
        try:
            newick = mrca_node.subtree_newick
        except Exception as x:
            errstream.write('Could not fetch the subtree. Error: {}\n'.format(str(x)))
        else:
            errstream.write('The newick representation of the subtree rooted at this node is:\n')
            output.write('{}\n'.format(newick))

induced_newick = tree_of_life.induced_subtree(ott_ids=id_list)['subtree']
'''