#!/bin/sh
test_path="${MODCELLHPC_PATH}/test/2"

problem_path="${MODCELLHPC_PATH}/cases/ecoli-core/"
prodnet_path="${MODCELL2_PATH}/problems/ecoli-core/prodnet.mat"
params_file="${test_path}/test2.params"
ini_pop_file=""
output_file="${test_path}/test2-out.pop"
output_file_csv="${test_path}/test2-out.csv"

# Run modcell
eval "${MODCELLHPC_PATH}/src/modcell $problem_path $params_file $output_file $ini_pop_file" || exit

# Convert ouput
eval "${MODCELLHPC_PATH}/io/pop2csv.py $problem_path $output_file -o $output_file_csv" || exit

# Check with matlab
temp_script=$(mktemp)
echo "cd ${test_path}" >> $temp_script
echo "test_objectives(\"${output_file_csv}\", \"${prodnet_path}\")" >> $temp_script
eval "${MATLAB_BIN} -nodesktop -nodisplay -sd ~/wrk/s/matlab < $temp_script"

