#!/bin/sh

cd ecoli-core

prodnet_path='~/wrk/s/matlab/modcell2/problems/ecoli-core/prodnet.mat'
temp_script=$(mktemp)

echo "cd  ~/wrk/s/matlab/modcell-hpc/cases/ecoli-core" >> $temp_script
echo "prodnet2txt(\"${prodnet_path}\")" >> $temp_script

~/progs/matlab/bin/matlab -nodesktop -nodisplay -sd ~/wrk/s/matlab < $temp_script

