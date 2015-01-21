#!/bin/bash

# move to analysis
cd ../analysis/

# run preprocessing of the data
R --vanilla < ../src/fat.preproc.R

# run im qtl analysis & permutations
R --vanilla < ../src/fat.qtl.R

# run cim qtl analysis
R --vanilla < ../src/fat.cim.R

# to run plots succesfully, you require the 3x3 rCIM factorial design
R --vanilla < ../src/fat.plots.R
