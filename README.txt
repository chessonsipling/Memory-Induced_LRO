This is the documentation for the manuscript "Memory-induced long-range order in dynamical systems" by C. Sipling and M. Di Ventra. Our work's thesis, background information, general analysis, and methods can be found in that paper (and/or in its included references).

Upon downloading all the *.py files in this repo, exe.py is ready to be run from the primary directory. This will create a new subdirectory in which spin lattice systems are simulated and avalanche distributions are extracted/plotted. Individual functions can be commented out depending on which specific features you'd like to be performed.

direc_copier.py can be run from the primary directory as well, which will create a new (filtered) subdirectory with all the figures and avalanche statistics from the original (raw) subdirectory, but without the storage-intensive raw data files (the explicit spin and memory dynamics). Make sure to change the names of the desired raw and filtered subdirectories in direc_copier.py itself.

See parameter_notes.txt for a description of all parameters which can be tuned, listed in the order they are defined in exe.py