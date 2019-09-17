#!/bin/sh

test_path="${MODCELLHPC_PATH}/test/io_1"
problem_path="${MODCELLHPC_PATH}/cases/ecoli-core/"
input_pop="${MODCELLHPC_PATH}/test/2/out.pop"
output_csv="${test_path}/pop2csv_out.csv"
output_pop="${test_path}/csv2pop_out.pop"

#eval "${MODCELLHPC_PATH}/io/pop2csv.py $problem_path $input_pop -o $output_csv"

eval "${MODCELLHPC_PATH}/io/csv2pop.py $problem_path $output_csv -o $output_pop"
