#! /usr/bin/env python

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
    refname = orig.readline()
    header = new.readline()
    out.write(header)
    # nuc_count = 0
    i = 0
    refnuc = 1
    snps.write("Differences between ref: {} and {}".format(refname, header))
    while refnuc:
        # nuc_count+=1
        refnuc = orig.read(1)
        # print(refnuc)
        # print(nuc_count)
        if refnuc == '-':
            print("Gap_found")
            out.write('-')
        else:
            if refnuc.lower() in set(['a', 'c', 'g', 't', 'n', 'r', 'y', 'm', 's', 'w', 'k', 'v', 'd', 'h', 'b']):
                newnuc = new.read(1)
                if newnuc and (newnuc != '\n'):
                    i+=1
                    if refnuc.lower() != newnuc.lower():
                        # print("position {}, ref {}, alt {}\n".format(i,refnuc,newnuc))
                        snps.write("position {}, ref {}, alt {}\n".format(i,refnuc,newnuc))
                    out.write(newnuc)
                else:
                    out.write('n')
            else:
                if refnuc:
                    print(refnuc)
                    assert refnuc =='\n'
                    out.write(refnuc)
    new.close()
    orig.close()
    snps.close()
    out.close()

if __name__ == '__main__':
    main()
