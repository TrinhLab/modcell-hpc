# ModCell-HPC
ModCell-HPC is an implementation of the ModCell2 multi-objective strain design method for optimized performance to deal with many-objective problems. Additionally, modcell-hpc:
- Does not depend on proprietary programs such as Matlab to run, which among other advantages allows you to run the program easily on any HPC resources.
- It has ~1000 lines of source code distributed across  5 `.c` files and the only external dependency is the glpk libary. This makes the program easy to understand and modify.
- Follows the unix philosophy thus modcell-hpc performs a well defined task and uses text files as input and output interfaces, this makes programming extensions to prepare the input or analyze the output very simple.

If you use any part of this software please cite:
~~~
        Sergio Garcia and Cong T. Trinh.
	In preparation
~~~

For other relevant studies check [these](link/to/bibfile) references in bib format.

This repository mainly contains the modcell-hpc tool with small examples. For the results and data analysis tools used in the study above visit the [modcell-hpc-study](https://github.com/TrinhLab/modcell-hpc-study) repository.

## Usage
This tool and instructions are for unix-like OSs (e.g., Linux, MacOS, BSD), if you want to follow along on Windows the simplest solution would be to run it on a Linux virtual machine.

### Input files
The input files correspond to LP problems for each production network (in `.mps` format) and information about candidates (IDs of reactions that can be deleted). These can be generated with the assistance of [ModCell's Matlab implementation](https://github.com/TrinhLab/ModCell2), but that is not a necessity.

### Output files
The output of the main method corresponds to a plain text files with information about each individual in the Population (design variables, design objectives, etc). This file can be converted into a table that only preserves Pareto optimal solutions and is useful for further analysis using the program `./io/pop2csv.py`. The resulting table can be analyzed with the help the small programs provided in `./tools`

### Running modcell-hpc
Run the `modcell` binary (either the released version or compile it your self as described below), the only runtime dependency is Open MPI (or any other MPI implementation). For necessary arguments and available options run `modcell --help`.

You can use scripts here or in [modcell-hpc-study](https://github.com/TrinhLab/modcell-hpc-study). Note that these scripts used predefine environment variables that correspond to paths in your system. So edit the file `paths` accordingly and add it to your shell by executing `source paths`. This needs to be done for every new shell, so instead you can add a line like this to your `~/.profile` or shellrc:
`[ -f "/path/to/modcell-hpc/paths" ] && source "/path/to/modcell-hpc/paths"`

## Project structure
- __src__ Contains modcell-hpc source code.
- __io__ Small programs that help generate input and process output of modcell-hpc.
- __tools__ Small programs to analyze the `.csv` output table containing genetic manipulations and associated design objectives for the identified Pareto optimal solutions.
- __cases__ Input files that include metabolic networks and candidate reactions.
- __test__  Programatic tests to ensure intended functionality.

## Filtering solutions
 The `.csv` file produced by `./io/pop2csv.py` noted above will remove duplicated, non-Pareto optimal, and solutions that violate alpha constraints. However, if module reactions are used, the output can still be improved further:

- __Minimizing module reactions__: The formulation of the modcell multi-objective optimization problem allows _module reactions_ that while satisfying the imposed constraints, do not help increase the objective for a given production network. These are denoted as futile module reactions and can be removed from an output population file by running `modcell` with the `--minimize_modules` flag. See `./test/min_1/test.sh` for an example.

- __Coalescing module reactions__: Solutions with the same reaction deletions but different module reactions can be combined into a superior solution. For further explanation and to perform this action check `./tools/mc_coalesce_modules`. See `./test/tools_1/test.sh` for an example.

## Compiling
You can use the provided Makefile. The following dependencies are needed:

- [GLPK](https://www.gnu.org/software/glpk/) (e.g., install in Arch Linux: `sudo pacman -S glpk`)
- MPI (Your favorite implementation, e.g., install in Arch Linux `sudo pacman -S openmpi`)

### Compilation examples
Run within `modcell-hpc/src`:

	- Compile with dynamic link (**Required** flag if you do not compile GLPK as explained below).
		- `make link=dynamic`
	- Compile with optimize flag:
		- `make flags=optimize`
	- Compile with optimize flag and static linking (will increase executable size, but remove the need for glpk installation, might also be a bit faster):
		- `make flags=optimize link=static`
	- Compile portable version (does not optimize to specific architecture and statically links glpk)
		- `make flags=portable`

### Compiling GLPK
To statically link glpk you probably need to compile it locally (unless your installation includes static libraries which is not common). A nice thing about this is that compilation flags can be tuned. The steps are as follows:

1. Obtain glpk:
	- `wget ftp://ftp.gnu.org/gnu/glpk/glpk-4.65.tar.gz`
	Verify download (optional):
		- `wget ftp://ftp.gnu.org/gnu/glpk/glpk-4.65.tar.gz.sig`
		- `gpg --recv-key 0xD17BF2305981E818 # Can be determined from running gpg glpk-4.65.tar.gz.sig , obviously you should check the key and make sure is still valid`
		- `gpg glpk-4.65.tar.gz.sig`

2. extract tarball and cd into glpk-X-YY:
	- `tar -xzvf glpk-4.65.tar.gz`
	- `cd glpk-4.65`
4. `./configure`
5. `make clean && make CFLAGS="-O3"`
6. The desired file is found under `glpk-X-YY/src/.libs/libglpk.a`

## Notes

### How does it work?
- The MOEA of choice is the proven NSGA-II.
- The ``flux balance analysis'' linear programming problems that determine metabolic fluxes are solved using GLPK.

### Why not use existing GA/MOEA libraries?
There are many libraries in various fast languages to do GA/MOEA. However, the very particular specifications and overall small size of this program have led me to conclude that avoiding these libraries will be beneficial for simplicity and optimization.

### Should I use the Matlab version or this version of modcell?
It depends on your preference. This version is much faster and requires less memory, even if you just run it on a personal computer. You can also simulate more cores than you have in your machine by simply adding the `--oversubscribe` flag to `mpiexec`, which can lead to better solution diversity.

### How do I run this program on a high performance computer (HPC)?
Using a HPC can sound daunting if you have never done it. However, it is essentially the same as running the program in your local machine, except that the hardware is much better. If you already have access to such system through your university, read their basic usage instructions, often MPI is installed and you will just need to write a small script to submit your job. Alternatively, you can access HPC through commercial providers, which often have free trial periods. Commercial providers offer high scalability and do not require waiting to run your jobs, however, they require additional setup steps. In this case, the simplest solution is to use one machine with many cores, or if you want a proper cluster you can try some automated solutions like [elasticluster](https://github.com/elasticluster/elasticluster).

### I ran into an error
The errors returned by C programs can be very uninformative and difficult to diagnose. The program does a good amount of consistency checking, but some random character in the arguments resulting in an obscure error message in the middle of it of the program. The best approach is to double check all input carefully before trying to diagnose a bug in the code itself. Of course if you are stuck or found an actual bug you should open an [issue](https://github.com/TrinhLab/modcell-hpc/issues/new).


### Credits
- [GLPK](https://www.gnu.org/software/glpk/) is used to solve LP problems.
- [PCG Random Number Generator](http://www.pcg-random.org/) is used to obtain quickly generated and uniformly distributed random numbers.
- [UTLIST and UTHASH](https://troydhanson.github.io/uthash/) are used for linked list and hash table data structures, respectively.
- [MPS format exporting tool](https://www.mathworks.com/matlabcentral/fileexchange/19618-mps-format-exporting-tool) is used to convert from the production networks generated by the [ModCell2 Matlab package](https://github.com/TrinhLab/ModCell2) to plain text .mps format.
