#! /usr/bin/env python3
import argparse
import os
import re
import subprocess as sp

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--test_all')
    return parser.parse_args()

def main():
    TEST_DIR = os.getcwd()
    EP_DIR = os.path.realpath(TEST_DIR + "/..")

    # Test the help menu of EP to assert it returns what we expect
    test_help_menu(TEST_DIR, EP_DIR)

    # Test a basic EP run
    # Tests the following flags
    # -a, -d, -1, -2, -u Align, -o
    test_basic_run(EP_DIR)




def test_help_menu(test_dir, ep_dir):

    ep_help_menu_expected_output = "Correctversionofbfctoolsfound.\nCorrectversionofsamtoolsfound.\nseqtkfound\nbwa-mem2found\nfastxtoolkitfound\nvcfutils.plfound\n\n\nExtensiphyisaprogramforquicklyaddinggenomicsequencedatatomultiplesequencealignmentsandphylogenies.\nViewtheREADMEformorespecificinformation.\nInputsaregenerallyamultiplesequencefileinfastaformatandadirectoryof\nFastqpaired-endreadsequences.\n\n\nEXAMPLECOMMAND:\n\n/path/to/extensiphy.sh-uALIGN-a/path/to/alignment_file-d/path/to/directory_of_reads[anyotheroptions]\n\nREQUIREDFLAGS\n(-a)alignmentinfastaformat,\n(-d)directoryofpairedendfastqreadfilesforallquerytaxa,\n(-u)produceonlyanupdatedalignmentorperformfullphylogeneticestimation(ALIGNorPHYLO)(DEFAULT:ALIGN)\n\n\nOPTIONALFLAGS\n(-t)treeinNewickformatproducedfromtheinputalignmentthatyouwishtoupdatewithnewsequencesorspecifyNONEtoperformnewinference(DEFAULT:NONE),\n(-m)alignmenttype(SINGLE_LOCUS_FILES,PARSNP_XMFAorCONCAT_MSA)(DEFAULT:CONCAT_MSA),\n(-o)directorynametoholdresults(DEFAULT:createsEP_output),\n(-i)cleanupintermediateoutputfilestosaveHDspace(Options:CLEAN,KEEP)(DEFAULT:KEEP),\n(-r)Selectedareferencesequencefromthealignmentfileforreadmappingorleaveasdefaultandarandomreferencewillbechosen(DEFAULT:RANDOM),\n(-p)numberoftaxatoprocessinparallel,\n(-c)numberofthreadspertaxonbeingprocessed,\n(-e)setread-typeassingleend(SE)orpair-end(PE)(DEFAULT:PE)\n(-1,-2)suffix(ex:R1.fastqorR2.fastq)forbothsetsofpairedendfiles(DEFAULTS:R1.fqandR2.fq),\n(-g)outputformat(CONCAT_MSAorSINGLE_LOCUS_FILES)(DEFAULT:CONCAT_MSA),\n(-s)specifythesuffix(.fa,.fasta,etc)(DEFAULT:.fasta),\n(-b)bootstrappingtreeONorOFF(DEFAULT:OFF)\n\n\nifusingsinglelocusMSAfilesasinput,\n(-f)csvfilenametokeeptrackofindividuallociwhenconcatenated(DEFAULT:loci_positions.csv),\n(-n)Setsizeoflocisizecutoffusedasinputoroutput(Options:intnumber)(DEFAULT:700)\n"

    no_lines_expected_output = ep_help_menu_expected_output.replace("\n", "")

    ep_help_menu = sp.Popen([ep_dir + '/extensiphy.sh', '-h'], stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = ep_help_menu.communicate()
    clean_output = str(stdout.decode("utf-8")).replace(" ", "").replace("\n","")

    assert no_lines_expected_output == clean_output
    print("EP help mentu test: PASSED")


def test_basic_run(ep_dir):

    output_location = ep_dir + "/tests/EP_output"

    ep_run = sp.Popen([ep_dir + "/extensiphy.sh", "-a", ep_dir + "/testdata/combo.fas", "-d", ep_dir + "/testdata", "-1", "_R1.fq", "-2", "_R2.fq", "-u", "ALIGN", "-o", output_location], stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = ep_run.communicate()

    # print(stdout.decode("utf-8"))
    find_alignment(output_location)

    ep_remove = sp.Popen(["rm", "-r", output_location], stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = ep_remove.communicate()


def find_alignment(DIR):
    alignment_file_types = ["aln", "fas", "fasta"]
    alignment_path = DIR + "/RESULTS/"
    files = os.listdir(alignment_path)

    for file in files:
        split_file_name = file.split(".")
        suffix_check = any(split_file_name[1] in x for x in alignment_file_types)
        if suffix_check:
            check_alignment(alignment_path, file)

def check_alignment(PATH, FILE):
    with open (PATH + "/" + FILE) as fasta:
        for line in fasta:
            print(line)

if __name__ == '__main__':
    main()
