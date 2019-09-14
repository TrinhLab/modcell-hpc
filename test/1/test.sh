#!/bin/sh
# Test dependent
TEST_N="1"
problem_path="${MODCELLHPC_PATH}/cases/ecoli-core/"
prodnet_path="${MODCELL2_PATH}/problems/ecoli-core/prodnet.mat"
ini_pop_file=""

#
test_path="${MODCELLHPC_PATH}/test/${TEST_N}"
params_file="${test_path}/Test.params"
ini_pop_file=""
output_file="${test_path}/out.pop"
output_file_csv="${test_path}/out.csv"
log_path="${test_path}/out.log"

date > "$log_path"
# Run modcell
eval "${MODCELLHPC_PATH}/src/modcell $problem_path $params_file $output_file $ini_pop_file" | tee -a "$log_path"|| exit

# Convert ouput
eval "${MODCELLHPC_PATH}/io/pop2csv.py $problem_path $output_file -o $output_file_csv" | tee -a "$log_path" || exit

# Check with matlab
temp_script=$(mktemp)
echo "cd ${test_path}" >> $temp_script
echo "test_objectives(\"${output_file_csv}\", \"${prodnet_path}\")" >> $temp_script
eval "${MATLAB_BIN} -nodesktop -nodisplay -sd ~/wrk/s/matlab < $temp_script" | tee -a "$log_path"

