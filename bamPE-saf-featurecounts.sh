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
##Count mapped reads for genomic features.

##Input $1: Output location.
##Input $2: Sorted and indexed paired-end BAM file (.bam).
##Input $3: Simplified Annotation Format file (.saf) for exons (with gene ID).
##Output: Counts (.txt) and counts summary (.summary).

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
##    -J bamPE-saf-featurecounts-${<INPUTFILE>##*/} \
##    bamPE-saf-featurecounts.sh <OUTPUT LOCATION> <INPUT BAM FILE> <INPUT SAF FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "bamPE-saf-featurecounts.sh") 

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
module load subread

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTBAMFILE=$(readlink -f $2)
INPUTSAFFILE=$(readlink -f $3)

INPUTBAMFILELOCATION=${INPUTBAMFILE%/*}
INPUTBAMFILENAME=${INPUTBAMFILE##*/}

############################################################################
##OUTPUT:

INPUTBAMFILEPREFIX=$(echo ${INPUTBAMFILENAME} | sed 's/\.bam$//' | sed 's/-job[0-9].*$//')
OUTPUTFILENAME=$(echo "${INPUTBAMFILEPREFIX}.featurecounts-job${JOBID}.txt") 
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}") 

############################################################################
##ACTIONS:

##Count features.
featureCounts -T 10 -F SAF -p --countReadPairs -a ${INPUTSAFFILE} -o ${OUTPUTFILE} ${INPUTBAMFILE}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTBAMFILE}
Input file: ${INPUTSAFFILE}
Output file: ${OUTPUTFILE}
Output file: ${OUTPUTFILE}.summary
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