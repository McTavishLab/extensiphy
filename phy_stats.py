#! /usr/bin/python3

import sys

nucleotides = ['A', 'C', 'G', 'T']
gaps = ['-']

gap_count=0
nucleotide_count=0
with open(sys.argv[1], 'r') as fasta:
    for i in fasta:
        if ">" not in i:
            for j in i:
                if j in nucleotides:
                    nucleotide_count+=1
                elif j in gaps:
                    gap_count+=1

percent_gap = (gap_count / nucleotide_count) * 100
percent_gap = float(str(percent_gap)[:5])

print("Percentage of gaps in the selected alignment")
print(percent_gap)
print("Number of gaps in the selected alignment")
print(gap_count)
print("number of nucleotides in the selected alignment")
print(nucleotide_count)
