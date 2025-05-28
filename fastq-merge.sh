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
##Merge FastQ files.

##Input $1: Output location.
##Input $2: FASTQ file (*.fastq.gz or *.fq.gz). 
##Input $3: FASTQ file (*.fastq.gz or *.fq.gz). 
##Output: Merged FASTQ file (*.fastq.gz or *.fq.gz).

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
##    -J fastq-merge-${<INPUTFILE>##*/} \
##    fastq-merge.sh <OUTPUT LOCATION> <INPUT FASTQ FILE 1> <INPUT FASTQ FILE 2> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "fastq-merge.sh") 

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



############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE1=$(readlink -f $2)
INPUTFILE2=$(readlink -f $3)

INPUTFILE1LOCATION=${INPUTFILE1%/*}
INPUTFILE2LOCATION=${INPUTFILE2%/*}
INPUTFILE1NAME=${INPUTFILE1##*/}
INPUTFILE2NAME=${INPUTFILE2##*/}

############################################################################
##OUTPUT:

INPUTFILE1PREFIX=$(echo ${INPUTFILE1NAME} | sed 's/\.fastq.gz$//' | sed 's/\.fq.gz$//' | sed 's/-job[0-9].*$//')
INPUTFILE2PREFIX=$(echo ${INPUTFILE2NAME} | sed 's/\.fastq.gz$//' | sed 's/\.fq.gz$//' | sed 's/-job[0-9].*$//')

##Create consensus file name.
CONSENSUSFILEPREFIX=""
while read C ; do
    C1=$(echo ${C} | cut -f1)
    C2=$(echo ${C} | cut -f2)
    if [[ "${C1}" == "${C2}" ]] ; then
        CONSENSUSFILEPREFIX+="${C1}"
    else
        CONSENSUSFILEPREFIX+="X"
    fi
done <<< "$(paste -d'\t' <(echo ${INPUTFILE1PREFIX} | grep -o .) <(echo ${INPUTFILE2PREFIX} | grep -o .))"

OUTPUTFILENAME=$(echo "${CONSENSUSFILEPREFIX}.mergefastq-job${JOBID}.fastq.gz")
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}") 

############################################################################
##ACTIONS:

cat ${INPUTFILE1} ${INPUTFILE2} > ${OUTPUTFILE}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE1}
Input file: ${INPUTFILE2}
Output file: ${OUTPUTFILE}
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