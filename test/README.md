Different types of tests for modcell-hpc

# How to run
Either `./run_all` to run all tests or `./run_test <test_directory_name>` to run a specifi one, e.g. `./run_test 1`  `./run_test io_1`

# Descriptions

## Objective tests

These are integration tests checking:
1.  Objective values computed by modcell-hpc can be reproduced independently (by Matlab ModCell2)
2. Known cases produce expected solutions

Tests:
- 1 : Basic test with beta = 0
- 2 : Basic test with beta > 0
- 3 : Test reading from an initial population
- 4 : Basic test with beta = 0 + MPI
- 5 : Basic test with beta > 0 + MPI
- 6 : MPI test random migration topology

## Other tests

Test io:
- io_1 : tests csv2pop
- io_2 : tests that pop2csv correctly identifies individuals that violate constraints, duplicates, and non-dominated individuals in a  .pop file

Test module minimizer:
- min_1: Tests the module minimization mode of modcell-hpc

Test tools
- tools_1: Tests mc_colesce_modules
- tools_2: Tests mc_dropbelowcomp
