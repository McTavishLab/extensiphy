# Extensiphy Testing Suite

Welcome to Extensiphy's Testing Suite. The tests here will help you check Extensiphy's functionality if you suspect the program isn't working correctly. You can run each test individually or you can run them all together (recommended). You can run these tests as many times as you like as the output directories produced by each test run will be removed prior to new tests being run.

## Quick Start

Run this command to test each major aspect of Extensiphy at once.

```
./comprehensive_test.sh
```

Thats it. Once the tests are complete, you should hopefully see this output:

```
test help menu: PASSED
test simple alignment update: PASSED
test update alignment with single-end reads: PASSED
test alignment update and phylo build: PASSED
test alignment update and phylo update: PASSED
test alignment update and bootstrap phylo update: PASSED
test alignment update and phylo update for specific reference: PASSED
test build and update alignment from single locus files: PASSED
test output single locus alignment files: PASSED
```

If any of these tests have failed and show a `FAILED` result, there is a problem with Extensiphy  
To get help with your installatio of Extensiphy, contact: jtoscanifield@ucmerced.edu
