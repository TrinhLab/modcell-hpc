#!/bin/sh

display_usage() {
	printf "Merge several populations\n"
	printf "\t- First argument: Solution directory where the different .pop files are found \n"
	printf "\t- Second argument (optional): name of the targets, e.g., if the targets are \"out.pop_*\", the paramtere should be given \"out\". Default is \"out\"."
}

# -------- Argument parsing --------

if [ "${1}" = "-h" ]; then
	display_usage && exit 1
fi

solution_dir="$1"
name=${2:-out}

# -------- Main program --------
target="${solution_dir}${name}.pop"
md_pattern="METADATA|population_size|alpha|beta"
touch temp
for fullfile in "${target}_"* ; do
	cat $fullfile >> temp
done
#cat "$target_"* | grep -Ev "#ENDFILE|${md_pattern}" | sed -e "\$a#ENDFILE"  > "$target"
grep -Ev "#ENDFILE|${md_pattern}" temp | sed -e "\$a#ENDFILE"  > "$target" && rm temp
grep -E $md_pattern "${target}_0" | cat - out.pop > temp && mv temp "$target"
