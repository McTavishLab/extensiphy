#! /usr/bin/python

import os
import argparse
import shutil
import subprocess


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', action='store_true', help='delete certain files')
    parser.add_argument('-m', action='store_true', help='copy certain files to a location')
    parser.add_argument('--taxa_dir')
    parser.add_argument('--new_dir')
    parser.add_argument('--max_num')
    return parser.parse_args()

def main():
    args = parse_args()

    if args.m == True:

        files = os.listdir(args.taxa_dir)
        os.chdir(args.taxa_dir)
        num = int(args.max_num)

        for i in files:
            if os.path.isfile(i):

                splitter = i.split('_')
                taxon_num = splitter[1]
                try:
                    tx = float(taxon_num)
                except ValueError:
                    print("2nd value in file name is not a number. Use renamer script")
                if tx <= num:
                    print(i)
                    subprocess.call(['cp', args.taxa_dir+i, args.new_dir + '/'])
                #stdout, stderr = process.communicate()


    elif args.d == True:

        files = os.listdir(args.taxa_dir)
        os.chdir(args.taxa_dir)
        num = int(args.max_num)

        for i in files:
            if os.path.isfile(i):
                splitter = i.split('_')
                taxon_num = splitter[1]
                try:
                    tx = float(taxon_num)
                except ValueError:
                    print("2nd value in file name is not a number. Use renamer script")

                if tx <= num:
                    print(i)
                    subprocess.Popen(['rm', args.taxa_dir+i])
                    #stdout, stderr = process.communicate()

if __name__ == '__main__':
    main()
