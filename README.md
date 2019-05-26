(In progress)

# ModCell-HPC
This is a free/libre version of the ModCell2 multi-objective strain design method for optimized performance and freedom of use (i.e., does not depend on proprietary programs such as Matlab to run, which among other advantages allows you to run the program on any HPC resources). If you use any part of this software please cite:
~~~
        Sergio Garcia, Cong T. Trinh,
        Multiobjective strain design: A framework for modular cell engineering,
        Metabolic Engineering,
        Volume 51,
        2019,
        Pages 110-120,
        ISSN 1096-7176,
        https://doi.org/10.1016/j.ymben.2018.09.003.
~~~

Another relevant study  is:
~~~
(moea-comparison)
~~~

## Input files
The input files correspond to LP problems for each production network (in .mps format) and information about candidates (IDs of reactions that can be deleted). These can be generated with the assistance of [ModCell's Matlab implementation](https://github.com/TrinhLab/ModCell2), but that is not a necessity.

## Output files
The output of the main method corresponds to a plain text files with information about each individual in the Population (design variables, design objectives, etc). This file can be converted into a table that only preserves Pareto optimal solutions and is useful for further analysis using the program `src/formatpopulation.py`. The resulting table can be analyzed with the help the small programs provided in [modcell-designs](https://github.com/TrinhLab/modcell-designs/tree/master/src).

## Compiling
You can use the provided Makefile. For optimal performance adjust compilation flags for the processor architecture of choice. The following dependencies are needed:
- [GLPK](https://www.gnu.org/software/glpk/) (In Arch linux install with `sudo pacman -S glpk`)

# Tests

# Notes

## How does it work?
- The MOEA of choice is the proven NSGA-II.
- The ``flux balance analysis" linear programming problems that determine metabolic fluxes are solved using GLPK.

## Why not use existing GA/MOEA libraries?
There are many libraries in various fast languages to do GA/MOEA. However, the very particular specifications and overall small size of this program have led me to conclude that avoiding these libraries will be beneficial for simplicity and optimization.

## Credits
- [GLPK](https://www.gnu.org/software/glpk/) is used to solve LP problems.
- [PCG Random Number Generator](http://www.pcg-random.org/) is used to obtain fastly generated and uniformly distributed random numbers.
- [UTLIST and UTHASH](https://troydhanson.github.io/uthash/) are used for linked list and hash table data structures, respectively.
