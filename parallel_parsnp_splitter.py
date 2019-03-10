#! /usr/bin/env python
## ~/usr/bin/env/python

import sys
import re

print(sys.argv)

infile = sys.argv[1]

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

with open(infile) as fp:
    for result in re.findall(regex2, fp.read()):
        print(result)
