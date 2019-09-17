Different types of tests for modcell-hpc

# How to run
Make sure modcell environment variables are loaded
sh 1/test1.sh

# Descriptions

## Objective tests

These are integration tests checking:
- Objective values computed by modcell-hpc can be reproduced independently (by matlab modcell2)
- Known cases produce expected solutions

Cases:
- test_test_objectives : Check if `test_objectives.m` works as intended.
- 1 : Basic test with beta = 0
- 2 : Basic test with beta > 0
- 3 : Test reading from an initial population
- 4 : Basic test with beta = 0 + MPI
- 5 : Basic test with beta > 0 + MPI
- 6 : MPI test random migration topology
- 7 : MPI test migration policies
- 8 : Read one initial population with MPI
- 9 : Read multiple initial populations with MPI

- io_1 : tests pop2csv and csv2pop
