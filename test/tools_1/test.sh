#!/bin/bash

test_path="${MODCELLHPC_PATH}/test/tools_1"
input_csv="${test_path}/in.csv"
output_csv="${test_path}/out.csv"

eval "${MODCELLHPC_PATH}/tools/mc_coalesce_modules $input_csv -o $output_csv"

# Assert expected output:
printf "Assert output--------------------------------\n"
printf "Expected solutions:\t 2\n"
rawlines=$(wc -l < "$output_csv")
lines=$((rawlines-1))
printf "Computed solutions:\t %d\n" "$lines"
printf "Expected checksum:\t %s\n" "5e2d22956f683689c3aa78fa03dfbd39e29f85828505b5bd3d5909074d427250"
# if equal means ouput is identical, otherwise output could still be correct but changed due to different row order
printf "Computed checksum:\t %s\n" $(sha256sum "$output_csv" | awk '{print $1}')
