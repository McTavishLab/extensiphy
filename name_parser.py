#! /usr/bin/python

import sys
import re

alignment = sys.argv[1]

data = open(alignment,'r')
readdata = data.read()

strp = readdata.strip()
splt = strp.split('>')

newnames_oldnames = {}
OTU_count = 0

for i in splt:
    count = 0
    OTU_count+=1
    j = i.split('\n')
    for segment in j:
        count+=1
        if count==1 and len(segment)>2:
            newnames_oldnames['OTU_' + str(OTU_count)] = segment

print(newnames_oldnames)
