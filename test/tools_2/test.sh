#!/bin/bash

test_path="${MODCELLHPC_PATH}/test/tools_2"
input_csv="${test_path}/in.csv"
output_csv="${test_path}/out.csv"

eval "${MODCELLHPC_PATH}/tools/mc_dropbelowcomp $input_csv -o $output_csv"
