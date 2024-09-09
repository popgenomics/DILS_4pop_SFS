#!/usr/bin/bash
## launch DILS for 1 populations
## the provided argument is for --configfile, expecting the yaml file
module load pypy/2.7-5.10.0
module load snakemake/7.7.0
module load r/4.1.1
module load python/3.9
binpath="/shared/home/croux/softwares/DILS_4pop/bin"
snakemake --snakefile ${binpath}/Snakefile -p -j 150 --cores 150 --configfile ${1} --cluster-config ${binpath}/cluster.json --latency-wait 30 --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time}"
