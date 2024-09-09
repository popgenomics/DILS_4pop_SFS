#!/usr/bin/bash
## launch DILS for 1 populations
## the provided argument is for --configfile, expecting the yaml file
module unload snakemake
module unload python
module unload pypy
module load pypy/7.3.11
module load snakemake/5.3.0
module load R/4.2.2
module load python/3.11.0
module load openjdk/17.0.3
binpath="/home/camille.roux/DILS_4pop/bin"
snakemake --snakefile ${binpath}/Snakefile -p -j 150 --cores 150 --configfile ${1} --cluster-config ${binpath}/cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --partition=BRYOFIT"
