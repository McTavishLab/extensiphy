#! /usr/bin/python

import sys
import re

alignment = sys.argv[1]

data = open(alignment,'r')
readdata = data.read()

strp = readdata.strip()
splt = strp.split('>')

clean_splt = splt[1:]

newnames_oldnames = {}
OTU_count = 0

for i in clean_splt:
    count = 0
    OTU_count+=1
    j = i.split('\n')
    for segment in j:
        count+=1
        if count==1:
            newnames_oldnames['OTU_' + str(OTU_count)] = segment

OTU_name_dict = "OTU_name_dict.txt"
new_alignment = open("renamed_alignment.fa", 'w')

for key, value in newnames_oldnames.iteritems():
    for item in clean_splt:
        if re.match(value, item):
            replaces = re.sub(value, key, item)
            new_alignment.write(">")
            new_alignment.write(replaces)
            new_alignment.write('\n')

with open(OTU_name_dict, 'w') as f:
    for key, value in newnames_oldnames.items():
        f.write('%s:%s\n' % (key, value))
