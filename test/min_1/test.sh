#!/bin/bash

test_path="${MODCELLHPC_PATH}/test/min_1"
problem_path="${MODCELLHPC_PATH}/cases/ecoli-core/"
input_csv="${test_path}/out.csv"
output_pop="${test_path}/out.pop"
minimized_pop="${test_path}/min_out.pop"
output_csv="${test_path}/min_out.csv"
#-----

# Generate .pop file
eval "${MODCELLHPC_PATH}/io/csv2pop.py $problem_path $input_csv -o $output_pop"

# Run module minimizer with approparite parameters
objective_type="wgcp"
alpha=5
beta=2
population_size=$(grep -c "#INDIVIDUAL" < "$output_pop")
eval "${MODCELLHPC_PATH}/src/modcell $problem_path $minimized_pop --initial_population=$output_pop --objective_type=$objective_type --alpha=$alpha --beta=$beta --population_size=$population_size --minimize_modules"

# Generate final csv file
eval "${MODCELLHPC_PATH}/io/pop2csv.py $problem_path $minimized_pop -o $output_csv -a $alpha"

# Assert expected output:
printf "Assert output--------------------------------\n"
printf "Expected solutions:\t 7\n"
rawlines=$(wc -l < "$output_csv")
lines=$((rawlines-1))
printf "Computed solutions:\t %d\n" "$lines"
printf "Expected checksum:\t %s\n" "6eb28b37edac7515e931ea83ef71e0ae236beabf204f0cc19d5041a84c971059"
# if equal means ouput is identical, otherwise output could still be correct but changed due to different row order
printf "Computed checksum:\t %s\n" $(sha256sum "$output_csv" | awk '{print $1}')
