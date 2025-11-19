# My scripts

- With the exception of *\*stdout.sh* scripts, all sbatch *\*.sh* scripts require the following environment variables:

	- **PATHTOMYSUBMITTEDSCRIPTS**: Path to directory where a copy of the submitted script will be saved.
	- **PATHTOMYSLURM**: Path to directory where slurm files are saved.
	- **MYSLURMFILE**: ${PATHTOMYSLURM}/slurm-%J.out

- Slurm output file should be set to **-o ${MYSLURMFILE}** during sbatch submission.

- *Rscript.sh* requires the following additional environment variable:

	- **PATHTOMYSCRIPTS**: Path to the directory where scripts are located.

- *Rscript-tmpscript.sh* requires the following additional environment variable:

	- **PATHTOMYTMPSCRIPTS**

- Some *\*stdout.sh* scripts require the following environment variables:

	- **PATHTOPROJTRASH**
	- **PATHTOPROJTMP**
	- **PATHTOPROJSAFE**
	- **PATHTOPROJEXPORT**
	- **PATHTOPROJTEST**
