#! /usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import dendropy

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree_file')
    return parser.parse_args()


def main():
    args = parse_args()

    mle = dendropy.Tree.get(path=args.tree_file, schema='newick')
    mle_len = mle.length()
    print("tree length is", mle_len)
    multifurcating = lambda x: True if len(x.child_nodes()) == 2 else False
    #for edge in mle.preorder_edge_iter():
    #for node in mle.preorder_node_iter():
        #print("node = ",node)
        #print(node.child_edges())
    for nd in mle.postorder_node_iter(multifurcating):
        print(nd.description(0))
        print(nd.edge_length)

    #for edge in mle.postorder_edge_iter():
    #    if edge.length is None:
    #        edge.length = 0
    #    else:
    #        edge.length = float(edge.length)/mle_len
    #print(mle.as_string(schema="newick"))
    #print(mle_len)


if __name__ == '__main__':
    main()
