#! usr/bin/env python

import sys

fq=sys.argv[1]
fa=sys.argv[2]
query=sys.argv[3]

fql=open(fq).readlines()


assert fql[0].startswith('@')

fal=open(fa,'w')
for lin in fql:
	if lin.startswith('+'):
		break
	if lin.startswith('@'):
		fal.write('>{}\n'.format(query))
	else:
		fal.write(lin)

fal.close


