#! /usr/bin/python3

import sys

degen = ['I', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D']
nucleotides = ['A', 'C', 'G', 'T']
gaps = ['-']
Ns = ['N']

n_count=0
degen_count=0
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
                elif j in degen:
                    degen_count+=1
                elif j in Ns:
                    n_count+=1


percent_gap = (gap_count / nucleotide_count) * 100
percent_gap = float(str(percent_gap)[:5])

print("Percentage of gaps in the selected alignment")
print(percent_gap)
print("Number of gaps in the selected alignment")
print(gap_count)
print("Number of nucleotides in the selected alignment")
print(nucleotide_count)
print('Number of degenerate nucleotides (non-Ns) in the selected alignment')
print(degen_count)
print('Number of Ns in selected alignment')
print(n_count)
