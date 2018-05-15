#!/usr/bin/env python

__DESCRIPTION__ = '''
Add indels to new consensus sequence to align to other taxa

align_consensus --gapped-ref file.fasta --consensus file.fas

outputs aligned_cns.fas
'''


import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gapped-ref', default='/dev/stdin')
    parser.add_argument('--consensus', default='/dev/stdout')
    parser.add_argument('--outfile', default="aligned_cns.fas")
    return parser.parse_args()


def main():
    args = parse_args()
    out = open(args.outfile, 'w')
    snps = open("snps.txt", 'w')
    orig = open(args.gapped_ref)
    new = open(args.consensus)
    orig.readline()
    header = new.readline()
    out.write(header)
    i = 0
    c = 1
    while c:
        c = orig.read(1)
        if c == '-':
            out.write('-')
        else:
            d = new.read(1)
            if d and (d != '\n'):
                i+=1
                if c.lower() != d.lower():
                    snps.write("poss {}, {}, {}".format(i,c,d))
                out.write(d)
            else:
                if c:
                    assert c == 'N' or c =='\n'
                    out.write(c)


if __name__ == '__main__':
    main()
