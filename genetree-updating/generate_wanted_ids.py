#!/usr/bin/env python
import sys
from dendropy import Tree
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

E_VALUE_THRESH = 0.04

Entrez.email = "ejmctavish@gmail.com"


study_id=sys.argv[1]
tree_id=sys.argv[2]
seqmat=sys.argv[3]
mattype=sys.argv[4]

ott_ncbi="../ott_ncbi"
study_id = "pg_2739"
tree_id = "tree6601"
seqmat = "tree6601.fas"
mattype = "fasta"
print study_id

#LSU ASC tree
ott_ncbi="../ott_ncbi"
study_id = "pg_873"
tree_id = "tree1679"
#seqmat = "ASC_LSU.nex"
seqmat = "ASC_LSU_nogap.fas"
#mattype = "nexus"
mattype="fasta"
print study_id

phy = Phylesystem()


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

#named_node = tree_of_life.mrca(ott_ids=[mrca_node.nearest_taxon.ott_id], wrap_response=True)
#newick = named_node.subtree_newick # NOTE: Excludes daughter taxa that  are considered non monophylo.


ott_to_ncbi = {}
ncbi_to_ott = {}
fi =open(ott_ncbi)
for lin in fi:
    lii= lin.split(",")
    ott_to_ncbi[int(lii[0])]=int(lii[1])
    ncbi_to_ott[int(lii[1])]=int(lii[0])


fi.close()

'''matches = re.finditer('ott(\d*)', newick) 
wanted_ncbi_ids = []
unfound = []
for item in matches:
    try:
        wanted_ncbi_ids.append(ott_to_ncbi[item.group(1)]
    except:
        unfound.append(item.group(1))
'''

#NOW: work through the alignment, blast each seq (EACH?)
#and make a list of genbacnk id's that match.

#OR ancestral states?
#blastn -remote -query fulltest3/cns.fa -db nr -max_target_seqs 1 -outfmt "6 qseqid staxids"
'''fi = open("/home/ejmctavish/ncbi/gi_taxid_nucl.dmp")
gi_taxid={}
for lin in fi:
    gi_taxid[lin.split()]
'''
 bsave= 

equery = "txid{}[orgn]".format(ott_to_ncbi[mrca_node.nearest_taxon.ott_id])

for i, record in enumerate(SeqIO.parse(seqmat, mattype)):
    result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"),  entrez_query=equery)
    save_file = open("asc_blast_{}.xml".format(i), "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()


gi_to_ncbi = {}
new_seqs={}

#for i, record in enumerate(SeqIO.parse(seqmat, mattype)):
for i in range(26):
    result_handle = open("asc_blast_{}.xml".format(i))
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                   new_seqs[int(alignment.title.split('|')[1])] = hsp.sbjct
                   gi_to_ncbi[int(alignment.title.split('|')[1])] = ''
##SET MINIMUM SEQUENCE LENGTH



 #This part is sort of crazy slooow
'''fi = open("/home/ejmctavish/ncbi/gi_taxid_nucl.dmp")
for lin in fi:
    if int(lin.split()[0]) in gi_to_ncbi.keys():
        lii=lin.split()
        gi_to_ncbi[int(lii[0])] = int(lii[1])

fi.close()'''

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


newseqs = "newseqs_asc.fasta"
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


subprocess.popen("papara", "-t", tree, "-s" alignment, "-q",  newseqs, "-n", "extended") 
                          #run RAXML EPA on the alignments

#placement
raxmlHPC -m GTRCAT -f v -s papara_alignment.multi_consensus -t ${WD}/$tree -n ${nam}_PLACE


#Full run with starting tree from placements
raxmlHPC -m GTRCAT -s papara_alignment.multi_consensus -t ${WD}/$tree -n update



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