#! /usr/bin/env python
## ~/usr/bin/env/python

import sys
import re

print(sys.argv)

infile = open(sys.argv[1], 'r')
header = []
i=0
for lin in infile:
    if lin.startswith('#'):
        header.append(lin)
        i+=1

infile.close()


infile = open(sys.argv[1], 'r')
spliton = infile.read()[i:].split('>') #Starts reading file after the header
infile.close()


taxoncount = header[1].strip().split(' ')[1]
locuscount = header[-1].strip().split(' ')[1]
taxon_key = {}

header = header[2:-1]
for count in range(0,len(header),4):
    idnum = header[count].split(' ')[1].strip()
    name = header[count + 1].split(' ')[1].lstrip('contigsworkingclean').rstrip('.fastq.fasta\n')
    taxon_key[idnum] = name

taxa = set()

seq_dict = {}
loci = {}

for seq_chunk in spliton:
    # print(seq_chunk)
    lines = seq_chunk.split('\n')
    taxid = lines[0].split(':')[0]
    try:
        int(taxid)
        if taxid != '0':
            taxa.add(taxid)
        m = re.search("cluster\d*", lines[0])
        locus = m.group()
        loci[locus] = []
        if taxid in seq_dict:
            seq_dict[taxid][locus] = lines[1:]
        else:
            seq_dict[taxid] = {}
            seq_dict[taxid][locus] = lines[1:]
    except:
        pass


assert set(seq_dict.keys()) == taxa
#assert len(loci) == int(locuscount)
assert len(taxa) == int(taxoncount)

## combiiiiine
locus_order = list(loci.keys())
locus_order.sort()

output = open("combo.fas","w")


for taxon in seq_dict:
    output.write(">ta{}\n".format(taxon_key[taxon]))
    seq = ""
    for locus in locus_order:
        locus_seq = "".join([i.strip("=") for i in seq_dict[taxon][locus]])
        loci[locus].append(len(seq))
        a = set(locus_seq)
        b = set(['a','g','c','t','A','C','G','T','-'])
        try:
            assert a.issubset(b)
        except:
            print(a)
            print(b)
        seq = seq+"{s}{n}".format(s=locus_seq,n="N"*100)
    output.write("{}\n".format(seq))

for locus in loci:
    assert(len(set(loci[locus]))==1)


output.close()

# adding this so git has something to commit...
