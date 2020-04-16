#! /usr/bin/env python
## ~/usr/bin/env/python

import sys
import re

print(sys.argv)

datafile = sys.argv[1]

infile = open(datafile, 'r')

# generate regex variables to find which taxon is associated with the number of this run
# and to count the total number of taxa
header = []
i=0
name_count = 0
seq_name_id = "##SequenceIndex " + sys.argv[2] + "\n"
seq_name_id_grabber = re.compile(seq_name_id)

taxa_count_regex = "##SequenceIndex " + "\d+"
taxa_count_compile = re.compile(taxa_count_regex)
taxa_counter = 0

# run regex, gather taxon name
for lin in infile:
    if lin.startswith('#'):
        # header.append(lin)
        count_taxa = re.search(taxa_count_compile, lin)
        grab_taxa_name = re.search(seq_name_id_grabber, lin)
        i+=1
        if count_taxa:
            taxa_counter+=1

        if grab_taxa_name:
            name_count = i

        if i == name_count + 1:
            clean_name = lin.replace("##SequenceFile ", "")
            header.append(clean_name)

infile.close()

# regex for start of taxon header
# this will be used to grab each chunk of a sequence
string_to_convert1 = ">" + sys.argv[2] + ":"
regex1 = re.compile(string_to_convert1)

# handles event that taxon being parsed is the last taxon in the alignment
if int(sys.argv[2]) == taxa_counter:
    ending_num = str(1)
    prime_regex = string_to_convert1 + "(.*?)" + "="
    regex2 = re.compile(prime_regex, re.S)
    #print(prime_regex)

# handles all taxa except the last taxon in the alignment
elif int(sys.argv[2]) != taxa_counter:
    ending_num = str(int(sys.argv[2]) + 1)
    prime_regex = string_to_convert1 + "(.*?)" + ">" + ending_num + ":"
    regex2 = re.compile(prime_regex, re.S)
    #print(prime_regex)

else:
    print("Problem with taxa counting")
    exit

# this is the ugly section that grabs the sequences, removes stuff
# and processes the file into a fasta file
seq = []
with open(datafile) as fp:
    for result in re.findall(regex2, fp.read()):
        seq.append(result)

almost_clean_seq = []

for item in seq:
    num_reg = "(^.*:p\d+\s)"
    reg_compile = re.compile(num_reg)
    string_search = re.findall(reg_compile, item)
    if string_search:

        new_seq = item.split(''.join(string_search))
        almost_clean_seq.append(new_seq[1])


newlines_left_seq = []
for seq in almost_clean_seq:
    new_seq = ''.join(seq).replace("\n", "")
    newlines_left_seq.append(new_seq)

clean_seq = []
for seq in newlines_left_seq:
    clean_seq.append(seq)
    clean_seq.append("N"*100)

clean_seq = ''.join(clean_seq)


# write everything to file for combining
new_seq_file = open("parsnp_chunk-" + sys.argv[2] + "-.fa",'w')
new_seq_file.write("\n")
new_seq_file.write(">")
new_seq_file.write(header[1])
new_seq_file.write(clean_seq)
new_seq_file.close()
