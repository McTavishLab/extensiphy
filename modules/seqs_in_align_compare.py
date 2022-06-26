#! /usr/bin/env python3
import os
import pandas as pd
import argparse
import multiprocessing

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--align_file', help='alignment you wish to investigate')
    # parser.add_argument('--cores', default=2, help='The number of cores you wish to use for this job')
    return parser.parse_args()


def main():
    args = parse_args()

    align_path = os.path.abspath(args.align_file)

    read_align = open(align_path, 'r').read()

    split_align = read_align.split('>')

    lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]


    taxon_list = []
    for num, item in enumerate(split_align):
        name_and_seq = item.split('\n', 1)
        if len(name_and_seq) > 1:
            taxon_list.append(name_and_seq[0])

    comparison_table = make_table(taxon_list)
    # print(comparison_table)

    segmented_list_of_lists = lol(split_align, 154)

    # print(segmented_list_of_lists[0])

    prc1 = multiprocessing.Process(target=compare_seqs, args=(segmented_list_of_lists[0], split_align, comparison_table, 1, ))
    prc2 = multiprocessing.Process(target=compare_seqs, args=(segmented_list_of_lists[1], split_align, comparison_table, 2, ))
    prc3 = multiprocessing.Process(target=compare_seqs, args=(segmented_list_of_lists[2], split_align, comparison_table, 3, ))
    prc4 = multiprocessing.Process(target=compare_seqs, args=(segmented_list_of_lists[3], split_align, comparison_table, 4, ))

    prc1.start()
    prc2.start()
    prc3.start()
    prc4.start()

    prc1.join()
    prc2.join()
    prc3.join()
    prc4.join()

    # print(prc1)
    # print(comparison_table)
    # dist_table = compare_seqs(split_align)

    # comparison_table.to_csv('tmp_table.csv')



def make_table(taxon_list):
    """
    Takes in a list of taxon names and makes an all against all table
    """
    df = pd.DataFrame(columns=taxon_list, index=taxon_list)

    return df


def compare_seqs(split_align_list, full_list, df, run_number):
    """
    Itterates over each sequence in the alignment and compares it to each other sequence in the alignment. \
    Records the number of differences between each sequence.
    """
    # output_dict = {}
    for num_1, item_1 in enumerate(split_align_list):
        print(num_1)
        name_and_seq_1 = item_1.split('\n', 1)
        # name_1 = name_and_seq_1[0]
        # seq_1 = name_and_seq_1[1]
        for num_2, item_2 in enumerate(full_list):
            name_and_seq_2 = item_2.split('\n', 1)
            # name_2 = name_and_seq_2[0]
            # seq_2 = name_and_seq_2[1]
            if len(name_and_seq_1) > 1 and len(name_and_seq_2) > 1:
                differing_bases = 0
                name_1 = name_and_seq_1[0]
                seq_1 = name_and_seq_1[1]
                name_2 = name_and_seq_2[0]
                seq_2 = name_and_seq_2[1]
                zipped_seqs = zip(seq_1, seq_2)
                for pair in zipped_seqs:
                    if pair[0] != pair[1]:
                        differing_bases+=1

                # output_dict[differing_bases] = [name_1, name_2]
                try:
                    df.at[name_1, name_2] = differing_bases

                except KeyError:
                    print("KEYERROR ", name_1 + ' ' + Name_2)

                df.to_csv('dist_table_' + str(run_number) + '.csv')

    return df















if __name__ == '__main__':
    main()
