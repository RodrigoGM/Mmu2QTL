#analysis/

This directory contains two shell scripts that runs QTL analyses for mouse body fat data contained in the data/ directory.

If you have a big enough computer, mainly lots of RAM, you can run these directly as

```
./fat.qtl.analysis.sh

./fat.rcim.analysis.sh
```

Otherswise, deploy with a cluster resource management system such as SLURM, lsf, and SGE.

For slurm :

```
sbatch --nodes=1 --ntasks-per-node=1 --mem-per-cpu=4000 ./fat.qtl.analysis.sh

sbatch --nodes=1 --ntasks-per-node=1 --mem-per-cpu=24000 --dependency afterok:<jobid> ./fat.rcim.analysis.sh
```


*System requirements* : for the analysis you will require R:Language and Environment, and the packages r/qtl, boa, coda, hdrcde, and RColorBrewer. You can install these as : 
```
install.packages(c("qtl", "boa", "coda", "hdrcde", "RColorBrewer"))
```
All the scripts run on a single core, and most should work with ± 4GB of RAM. However, scripts for the repicated CIM analyses will require more RAM ± 16 GB or more.

R/qtl allows for mutlithreading with the n.cluster argument (requires packages snow and rlecuyer). In addition lapply and sapply loops can be run in paralell with ```snow``` using ```parLapply``` or ```parSapply```. , For compatibility with most systems these options were not implemented.
