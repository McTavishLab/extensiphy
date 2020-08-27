#! /bin/bash

taxa_count=$(grep -oP ">(\d+):" parsnp.xmfa | tail -1 | sed -e 's/>//' | sed -e 's/://')

num_seq=$(seq -s ' ' 1 $taxa_count)

for i in $num_seq; do
	echo $i
done


