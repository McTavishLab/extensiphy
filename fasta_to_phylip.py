#! /usr/bin/env python

__DESCRIPTION__ = '''
Convert a FASTA alignment to Phylip format.

Dependenies: BioPython

fasta_to_phylip --input-fasta file.fasta --output-phy file.phy

'''

import dendropy
import argparse


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input-fasta', default='/dev/stdin')
    parser.add_argument('--output-phy', default='/dev/stdout')

    return parser.parse_args()

def main():

    args = parse_args()

    aln = dendropy.DnaCharacterMatrix.get(
        path=args.input_fasta,
        schema="fasta")

    aln.write(
        path=args.output_phy,
        schema="phylip",
        strict=False)
    
if __name__ == '__main__':
    main()
