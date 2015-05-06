#!/usr/bin/python
import sys

fastafi=sys.argv[1]
outfi=sys.argv[2]

ofi=open(outfi,'w')
fas=open(fastafi).readlines()

ofi.write(fas[0])

for i, char in enumerate(''.join([lin.strip() for lin in fas[1:]])):
	if i%80==0 and i>0:
		ofi.write('\n')
	ofi.write(char)

ofi.close()