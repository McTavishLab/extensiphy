README for function-compare-test

This test will run test reads through gon_phyling.sh (traditional phylogenetic analysis pipeline) and phycorder (rapid updating method)
and compare the trees output by both methods by providing the RF distance between them.

This test will also ensure that the trees produced by both pipelines have the topology that they should.

To Use:

make sure method_compare.sh, reset_compare_test.sh and tree_compare_results.py are executeable on your machine.

run: ./method_compare.sh to run the test

reset test space: ./reset_compare_test.sh 
