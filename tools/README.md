ModCell data analysis tools following the unix philosophy.
These are command line tools that perform specific actions on the Pareto front expressed in .csv format.

## Requirements
- Most scripts depend on general data analysis libraries such as `pandas`.
- `setcover.py` depends on `pyomo` (`pip install pyomo`) and an optimization solver (usually this will be installed at the system level, e.g. in Arch Linux `pacman -S glpk`).
