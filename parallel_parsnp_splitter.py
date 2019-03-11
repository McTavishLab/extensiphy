#! /usr/bin/env python
## ~/usr/bin/env/python

import sys
import re

print(sys.argv)

datafile = sys.argv[1]

infile = open(sys.argv[1], 'r')
header = []
i=0
for lin in infile:
    if lin.startswith('#'):
        header.append(lin)
        i+=1

infile.close()

string_to_convert1 = ">" + sys.argv[2] + ":"
regex1 = re.compile(string_to_convert1)

ending_num = str(int(sys.argv[2]) + 1)
prime_regex = string_to_convert1 + "(.*?)" + ">" + ending_num + ":"
regex2 = re.compile(prime_regex, re.S)
print(prime_regex)


seq = []
with open(datafile) as fp:
    for result in re.findall(regex2, fp.read()):
        seq.append(result)

almost_clean_seq = []

for item in seq:
    num_reg = ":p\d+"
    reg_compile = re.compile(num_reg)
    string_search = re.findall(reg_compile, item)
    if string_search:
        # print(string_search)
        new_seq = item.split(''.join(string_search))
        almost_clean_seq.append(new_seq[1])

newlines_left_seq = []
for seq in almost_clean_seq:
    new_seq = ''.join(seq).strip()
    newlines_left_seq.append(new_seq)

print(newlines_left_seq)
# clean_seq = newlines_left_seq.strip("\n")
# print(clean_seq)
