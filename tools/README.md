ModCell data analysis tools following the unix philosophy.
These are command line tools that perform specific actions on the Pareto front expressed in .csv format.
Notation:
- mc is appended to make this tools easier to call from your shell by typing mc and using tab autocomplete
- mc_dir_* are tools that act on a directory containing different pareto fronts.
- mc_help_* provide generate files useful for other tools.

## Requirements
- Most scripts depend on popular data analysis libraries such as `pandas` and `seaborn`.
- `mc_setcover` depends on `pyomo` (`pip install pyomo`) and an optimization solver (usually this will be installed at the system level, e.g. in Arch Linux `pacman -S glpk`).
- I don't usually run these on a virtualenv but I have included `requirements.txt` since I have experienced cobratoolbox breaking compatibility.

