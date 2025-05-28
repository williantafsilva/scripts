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
##Merge BAM files (.bam).

##Input $1: Output location.
##Input $2: BAM file (.bam).
##Input $3: BAM file (.bam).
##Output: Merged BAM file (.bam) and BAM index file.

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
##    -J bam-merge-${<INPUTFILE>##*/} \
##    bam-merge.sh <OUTPUT LOCATION> <INPUT BAM FILE 1> <INPUT BAM FILE 2>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "bam-merge.sh") 

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
INPUTFILE1=$(readlink -f $2)
INPUTFILE2=$(readlink -f $3)

INPUTFILE1LOCATION=${INPUTFILE1%/*}
INPUTFILE2LOCATION=${INPUTFILE2%/*}
INPUTFILE1NAME=${INPUTFILE1##*/}
INPUTFILE2NAME=${INPUTFILE2##*/}

############################################################################
##OUTPUT:

INPUTFILE1PREFIX=$(echo ${INPUTFILE1NAME} | sed 's/\.bam.*$//' | sed 's/-job[0-9].*$//')
INPUTFILE2PREFIX=$(echo ${INPUTFILE2NAME} | sed 's/\.bam.*$//' | sed 's/-job[0-9].*$//')

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

OUTPUTFILENAME=$(echo "${CONSENSUSFILEPREFIX}.mergebam-job${JOBID}.bam")
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}") 

############################################################################
##ACTIONS:

##Merge BAM files.
samtools merge -o ${OUTPUTFILE} ${INPUTFILE1} ${INPUTFILE2}

##Sort merged BAM file.
samtools sort -o ${OUTPUTFILE} ${OUTPUTFILE}

##Index merged BAM file.
samtools index -b ${OUTPUTFILE}

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