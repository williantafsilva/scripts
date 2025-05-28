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
##Sort and index BAM file.

##Input $1: Output location.
##Input $2: BAM file (.bam).
##Output: Sorted BAM file (.bam) and index file (.bai).

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
##    -J bam-sort-index-${<INPUTFILE>##*/} \
##    bam-sort-index.sh <OUTPUT LOCATION> <INPUT BAM FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "bam-sort-index.sh") 

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
module load samtools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

INPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.bam.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILENAME=$(echo "${INPUTFILEPREFIX}.sortbam-job${JOBID}.bam")
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}") 

############################################################################
##ACTIONS:

##Sort BAM file.
samtools sort -o ${OUTPUTFILE} ${INPUTFILE}

##Index BAM file.
samtools index -b ${OUTPUTFILE}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Output file: ${OUTPUTFILE}
Output file: ${OUTPUTFILE}.bai
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