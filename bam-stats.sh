#!/bin/bash -l
##Add -l to the shebang to inherit bash profile variables and configuration.
##Required environment variables: ${PATHTOMYSUBMITTEDSCRIPTS}, ${PATHTOMYSLURM}, ${PATHTOPROJTMP}, ${PROJECT_ID}, ${MYSLURMFILE}, ${MYEMAIL}.
############################################################################
################################## SCRIPT ##################################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Calculate genome-wide statistics and coverage of a BAM file.

##Input $1: Output location.
##Input $2: Sorted and indexed BAM file aligned to coordinates (.bam).
##Output: Genome-wide statistics and coverage data (.txt).

##Usage: 
##sbatch \
##    -A ${PROJECT_ID} \
##    -o ${MYSLURMFILE} \
##    -p shared \
##    -N 1 \
##    -n 10 \
##    -t 05-00:00:00 \
##    --mail-type=ALL \
##    --mail-user=${MYEMAIL} \
##    -J bam-stats-${<INPUTFILE>##*/} \
##    bam-stats.sh <OUTPUT LOCATION> <INPUT BAM FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "bam-stats.sh") 

############################################################################
##JOB ID:

RUNDATE=$(date +"%Y%m%d%H%M%S")
if [[ -z "${SLURM_JOB_ID}" ]] ; then JOBID=${RUNDATE} ; else JOBID=${SLURM_JOB_ID} ; fi 

############################################################################
##SUBMITTED SCRIPT COPY:

cat $0 > $(echo "${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}") 

############################################################################
##DIAGNOSTICS:

echo "$(date +"%Y-%m-%d @ %H:%M:%S"): 
Submission: ${SCRIPTNAME} $@
User (\$USER): ${USER}
Host name (hostname -f): $(hostname -f)
Operating system: $(uname -a)
PATH: ${PATH}
Job ID (\$SLURM_JOB_ID): ${SLURM_JOB_ID}"

############################################################################
##LOAD TOOLS:

module load bioinfo-tools
module load bamtools
module load samtools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

INPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.bam$//' | sed 's/-job[0-9].*$//')
OUTPUTFILE1NAME=$(echo "${INPUTFILEPREFIX}.bamstatsbamtools-job${JOBID}.txt")
OUTPUTFILE2NAME=$(echo "${INPUTFILEPREFIX}.bamstatssamtools-job${JOBID}.txt")
OUTPUTFILE3NAME=$(echo "${INPUTFILEPREFIX}.bamstatsindex-job${JOBID}.txt")
OUTPUTFILE4NAME=$(echo "${INPUTFILEPREFIX}.bamstatsgenomecoverage-job${JOBID}.txt")
OUTPUTFILE5NAME=$(echo "${INPUTFILEPREFIX}.bamstatsgenomecoverageplots-job${JOBID}.txt")
OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE1NAME}") 
OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE2NAME}") 
OUTPUTFILE3=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE3NAME}") 
OUTPUTFILE4=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE4NAME}") 
OUTPUTFILE5=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE5NAME}") 

############################################################################
##ACTIONS:

##General statistics.
bamtools stats -in ${INPUTFILE} > ${OUTPUTFILE1}
samtools stats ${INPUTFILE} > ${OUTPUTFILE2}
samtools idxstats ${INPUTFILE} > ${OUTPUTFILE3}
samtools coverage ${INPUTFILE} >  ${OUTPUTFILE4}

##Genome-wide coverage histogram.
samtools coverage --ascii --plot-depth --n-bins 25 ${INPUTFILE} > ${OUTPUTFILE5}
echo "
Average genome-wide coverage (samtools depth): " >> ${OUTPUTFILE5}
samtools depth -a ${INPUTFILE} | awk '{sum+=$3; sumsq+=$3*$3} END {print "Average coverage = ",sum/NR; print "Std. Dev. = ",sqrt(sumsq/NR - (sum/NR)**2)}' >> ${OUTPUTFILE5}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE5}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Output file: ${OUTPUTFILE1}
Output file: ${OUTPUTFILE2}
Output file: ${OUTPUTFILE3}
Output file: ${OUTPUTFILE4}
Output file: ${OUTPUTFILE5}
" >> $(echo "${OUTPUTLOCATION}/README.txt") 

	##COPY OF SCRIPT:
	cat $0 > $(echo "${OUTPUTLOCATION}/job${JOBID}.sh")
	##COPY OF SLURM FILE:
	cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTLOCATION}/job${JOBID}.out") 

else

	##COPY OF SCRIPT:
	cat $0 > $(echo "${OUTPUTLOCATION}/job${JOBID}.sh")
	##COPY OF SLURM FILE:
	cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTLOCATION}/job${JOBID}.err") 

fi

exit 0