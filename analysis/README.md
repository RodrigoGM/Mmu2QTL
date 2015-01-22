#analysis/

This directory contains two shell scripts that runs QTL analyses for mouse body fat data contained in the data/ directory.

If you have a big enough computer, mainly lots of RAM, you can run these directly as

```./fat.qtl.analysis.sh```

```./fat.rcim.analysis.sh```

Otherswise, deploy with a cluster resource management system such as SLURM, lsf, and SGE.

For slurm :

 ```sbatch --nodes=1 --ntasks-per-node=1 --mem-per-cpu=4000 ./fat.qtl.analysis.sh```
 
```sbatch --nodes=1 --ntasks-per-node=1 --mem-per-cpu=24000 ./fat.rcim.analysis.sh```

