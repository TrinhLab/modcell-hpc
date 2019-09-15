#!/bin/sh
# Test dependent
TEST_N="6"
problem_path="${MODCELLHPC_PATH}/cases/ecoli-core/"
prodnet_path="${MODCELL2_PATH}/problems/ecoli-core/prodnet.mat"
ini_pop_file=""

# Parameters
objective_type="wgcp"
alpha=5
beta=0
population_size=100
n_generations=100
seed=0
crossover_probability=0.8
mutation_probability=0.05
max_run_time=7200
migration_topology=1

#
test_path="${MODCELLHPC_PATH}/test/${TEST_N}"
output_file="${test_path}/out.pop"
output_file_csv="${test_path}/out.csv"


# Run modcell
eval "mpiexec -n 4 ${MODCELLHPC_PATH}/src/modcell $problem_path $output_file --initial_population=$ini_pop_file --objective_type=$objective_type --alpha=$alpha --beta=$beta --population_size=$population_size --n_generations=$n_generations --seed=$seed --crossover_probability=$crossover_probability --mutation_probability=$mutation_probability --max_run_time=$max_run_time --migration_topology=$migration_topology" || exit

# Convert ouput
eval "${MODCELLHPC_PATH}/io/popmerge.sh $test_path/" || exit

# Convert ouput
eval "${MODCELLHPC_PATH}/io/pop2csv.py $problem_path $output_file -o $output_file_csv" || exit

# Check with matlab
temp_script=$(mktemp)
echo "cd ${test_path}" >> $temp_script
echo "test_objectives(\"${output_file_csv}\", \"${prodnet_path}\")" >> $temp_script
eval "${MATLAB_BIN} -nodesktop -nodisplay -sd ~/wrk/s/matlab < $temp_script"

