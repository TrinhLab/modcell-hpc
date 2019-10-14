Different types of tests for modcell-hpc

# How to run
Make sure modcell environment variables are loaded
sh 1/test1.sh

# Descriptions

## Objective tests

These are integration tests checking:
- Objective values computed by modcell-hpc can be reproduced independently (by matlab modcell2)
- Known cases produce expected solutions

Integration test:
- 1 : Basic test with beta = 0
- 2 : Basic test with beta > 0
- 3 : Test reading from an initial population
- 4 : Basic test with beta = 0 + MPI
- 5 : Basic test with beta > 0 + MPI
- 6 : MPI test random migration topology

Test tests:
- test_test_objectives : Check if `test_objectives.m` works as intended.

Test io:
- io_1 : tests csv2pop
- io_2 : tests that pop2csv correctly identifies individuals that violate constraints, duplicates, and non-dominated individuals in a  .pop file

Test module minimizer:
- min_1 :

# Known issues
- `test_objectives.m` seems to be working incorrectly
