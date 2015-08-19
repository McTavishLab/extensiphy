#!/usr/bin/env python
import sys
from dendropy import Tree, DnaCharacterMatrix
from peyotl import gen_otu_dict, iter_node
from peyotl.manip import iter_trees
from peyotl.phylesystem.phylesystem_umbrella import Phylesystem
from peyotl.sugar import tree_of_life, taxonomy
from peyotl.nexson_syntax import extract_tree,  PhyloSchema
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Entrez
import re
import subprocess
import time

'''
#LSU ASC tree
ott_ncbi="../ott_ncbi"
study_id = "pg_873"
tree_id = "tree1679"
#seqmat = "ASC_LSU.nex"
seqaln = "mat1.fas"
#mattype = "nexus"
runname="asc_test"
mattype="fasta"
print study_id
'''

study_id=sys.argv[1]
tree_id=sys.argv[2]
seqaln=sys.argv[3]
mattype=sys.argv[4]
runname=sys.argv[5]


#Fixed values
E_VALUE_THRESH = 0.04
ott_ncbi="../ott_ncbi" #TODO config file
Entrez.email = "ejmctavish@gmail.com"



phy = Phylesystem()
n = phy.return_study(study_id)[0]
api_wrapper.study.get(study_id,tree=tree_id)

##This is a weird way to get the ingroup node, but I need the OTT ids anyhow.
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

newick = extract_tree(n, tree_id, PhyloSchema('newick', output_nexml2json = '1.2.1', content="tree", tip_label="ot:originalLabel"))
newick = newick.replace(" ", "_") #UGH

d = DnaCharacterMatrix.get(path=seqaln,
                           schema=mattype)
# make the taxon_namespace immutable, so the tree does not add
#   new labels...
d.taxon_namespace.is_mutable = True
tre = Tree.get(data=newick,
                schema="newick",
                taxon_namespace=d.taxon_namespace)

# get all of the taxa associated with tips of the tree, and make sure that
#   they include all of the members of the data's taxon_namespace...
treed_taxa = [i.taxon for i in tre.leaf_nodes()]
if len(treed_taxa) != len(d.taxon_namespace):
    missing = [i.label for i in d.taxon_namespace if i not in treed_taxa]
    emf = 'Some of the taxa in the alignment are not in the tree. Missing "{}"\n'
    em = emf.format('", "'.join(missing))
    raise ValueError(em)

elem_ord = n['nexml']['^ot:otusElementOrder'][0]

map_dict = {}
idsset =set()
for otu in n['nexml']['otusById'][elem_ord]['otuById']:
    orig = n['nexml']['otusById'][elem_ord]['otuById'][otu]['^ot:originalLabel']
    try:
        ott_id = n['nexml']['otusById'][elem_ord]['otuById'][otu]['^ot:ottId']
        if ott_id in idsset:
            ott_id=orig.replace(" ","_")
        idsset.add(ott_id)
    except:
        ott_id=orig.replace(" ","_")
    map_dict[orig.replace(" ","_")]=str(ott_id)

for taxon in d.taxon_namespace:
    taxon.label = map_dict[taxon.label]


tre.resolve_polytomies()
tre.write(path = "{}_random_resolve.tre".format(runname), schema = "newick", unquoted_underscores=True, suppress_edge_lengths=True)

d.write(path="{}_aln_ott.phy".format(runname), schema="phylip")
d.write(path="{}_aln_ott.fas".format(runname), schema="fasta")

ott_to_ncbi = {}
ncbi_to_ott = {}
fi =open(ott_ncbi)
for lin in fi:
    lii= lin.split(",")
    ott_to_ncbi[int(lii[0])]=int(lii[1])
    ncbi_to_ott[int(lii[1])]=int(lii[0])


fi.close()


equery = "txid{}[orgn]".format(ott_to_ncbi[mrca_node.nearest_taxon.ott_id])

for i, record in enumerate(SeqIO.parse(seqmat, mattype)):
    result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"),  entrez_query=equery)
    save_file = open("{}_{}.xml".format(runname,i), "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()



gi_to_ncbi = {}
new_seqs={}

for i, record in enumerate(SeqIO.parse(seqmat, mattype)):
    result_handle = open("{}_{}.xml".format(runname,i))
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                   new_seqs[int(alignment.title.split('|')[1])] = hsp.sbjct
                   gi_to_ncbi[int(alignment.title.split('|')[1])] = ''
##SET MINIMUM SEQUENCE LENGTH


for gi in gi_to_ncbi.keys():
    if gi_to_ncbi[gi] == '':
        try:
            handle = Entrez.esummary(db="nucleotide", id=gi, rettype="gb", retmode="text")
            record = Entrez.read(handle)
            handle.close()
            gi_to_ncbi[gi] = record[0]['TaxId']
            time.sleep(1)
            sys.stdout.write(".")
        except:
            sys.stdout.write("x")
            pass


newseqs = "{}.fasta".format(runname)
fi = open(newseqs,'w')
for gi in gi_to_ncbi:
        try:
            ott_id = ncbi_to_ott[gi_to_ncbi[gi]]
        except:
            print("ncbi taxon ID {} not in ott".format(gi_to_ncbi[gi]))
            continue
        if ott_id not in ottids: # only adding seqs we don't have taxa for
                ottids.append(ott_id)
                fi.write(">{}\n".format(ncbi_to_ott[gi_to_ncbi[gi]]))
                fi.write("{}\n".format(new_seqs[gi]))
                print("success {}".format(ott_id))
        if ott_id in ottids: 
                print("ncbi taxon ID {} already in tree".format(gi_to_ncbi[gi]))


fi.close()

#Now get all descendants of that node in the synth tree. Or the taxonomy? better closest tax node above.
#NOW - > Need the real alignement
#aln = 
#and a RESOLVED TREE?! w/ BL?!


subprocess.Popen(["papara", "-t","{}_random_resolve.tre".format(runname), "-s", "{}_aln_ott.phy".format(runname), "-q",  newseqs, "-n", "extended"]) 
                          #run RAXML EPA on the alignments

#placement
subprocess.Popen(["raxmlHPC", "-m", "GTRCAT", "-f", "v", "-s", "papara_alignment.extended", "-t","{}_random_resolve.tre".format(runname), "-n", "{}_PLACE".format(runname)])


placetre = Tree.get(path="RAxML_labelledTree.{}_PLACE".format(runname),
                schema="newick")

for taxon in placetre.taxon_namespace:
    if taxon.label.startswith("QUERY"):
        taxon.label=taxon.label.split()[1]
    taxon.label=taxon.label.replace(" ", "_") #Look for real dendropy way!!


#NEEDS ASSERT THAT NAME SPACES ARE SAME!

placetre.resolve_polytomies()
placetre.write(path = "{}_place_resolve.tre".format(runname), schema = "newick", unquoted_underscores=True)
#Full run with starting tree from placements

subprocess.Popen(["raxmlHPC", "-m", "GTRCAT", "-s", "papara_alignment.extended", "-t","{}_place_resolve.tre".format(runname), "-p", "1", "-n", "{}".format(runname)])
