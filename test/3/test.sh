#!/bin/sh
# Test dependent
TEST_N="3"
problem_path="${MODCELLHPC_PATH}/cases/ecoli-core/"
prodnet_path="${MODCELL2_PATH}/problems/ecoli-core/prodnet.mat"
ini_pop_file=""

#
test_path="${MODCELLHPC_PATH}/test/${TEST_N}"
params_file="${test_path}/Test.params"
ini_pop_file=""
output_file="${test_path}/out.pop"
output_file2="${test_path}/out2.pop"
output_file_csv="${test_path}/out2.csv"

# Run modcell first
eval "${MODCELLHPC_PATH}/src/modcell $problem_path $params_file $output_file $ini_pop_file" || exit

# Run modcell again with initial population
eval "${MODCELLHPC_PATH}/src/modcell $problem_path $params_file $output_file2 $output_file" || exit

# Convert ouput
eval "${MODCELLHPC_PATH}/io/pop2csv.py $problem_path $output_file2 -o $output_file_csv" || exit

# Check with matlab
temp_script=$(mktemp)
echo "cd ${test_path}" >> $temp_script
echo "test_objectives(\"${output_file_csv}\", \"${prodnet_path}\")" >> $temp_script
eval "${MATLAB_BIN} -nodesktop -nodisplay -sd ~/wrk/s/matlab < $temp_script"

