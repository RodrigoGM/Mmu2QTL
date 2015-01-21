#!/bin/bash

# run preprocessing of the data
R --vanilla < ../src/fat.rcim.preproc.R

# run rcim in all step sizes
# each script loops through the three window sizes
R --vanilla < ../src/fat.rcim_s0.R
R --vanilla < ../src/fat.rcim_s1.R
R --vanilla < ../src/fat.rcim_s05.R

# summarizing the individual runs
R --vanilla < ../src/fat.rcim_smry_s0,1,05.R

# to run plots succesfully, you require the 3x3 rCIM factorial design
R --vanilla < ../src/fat.plots.R
