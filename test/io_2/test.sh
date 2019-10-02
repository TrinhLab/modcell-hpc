#!/bin/sh

test_path="${MODCELLHPC_PATH}/test/io_2"
problem_path="${MODCELLHPC_PATH}/cases/ecoli-core/"
input_pop="${MODCELLHPC_PATH}/test/io_2/crafted.pop"

eval "${MODCELLHPC_PATH}/io/pop2csv.py $problem_path $input_pop  -a 5"

echo "Expected: 5 solutions: 1 solution violates alpha, 2 are duplicate, and out of the 3 remaining, 2 have same objectives, and the 1 is dominated"

