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
##Input $2: List (.txt) with two tab-separated columns of FASTQ files to merge. 
##Output: Merged FASTQ files (*.fastq.gz or *.fq.gz).

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
##    -J list-fastq-merge-${<INPUTFILE>##*/} \
##    list-fastq-merge.sh <OUTPUT LOCATION> <INPUT FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "list-fastq-merge.sh") 

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
INPUTLIST=$(readlink -f $2)

INPUTLISTLOCATION=${INPUTLIST%/*}
INPUTLISTNAME=${INPUTLIST##*/}

############################################################################
##OUTPUT:

OUTPUTLIST=$(echo "${OUTPUTLOCATION}/files.mergefastq-job${JOBID}.txt")
ERRORFILE=$(echo "${OUTPUTLOCATION}/error.mergefastq-job${JOBID}.txt")

############################################################################
##ACTIONS:

echo -e "INPUTFILE1\tINPUTFILE2\tOUTPUTFILE\tERROR?" > ${OUTPUTLIST}
cat ${INPUTLIST} | while read LINE ; do

	#######################################
	##INPUT:
	
	INPUTFILE1=$(echo ${LINE} | cut -f1)
	INPUTFILE2=$(echo ${LINE} | cut -f2)

	#######################################
	##OUTPUT:

	INPUTFILE1PREFIX=$(echo ${INPUTFILE1##*/} | sed 's/\.fastq.gz$//' | sed 's/\.fq.gz$//' | sed 's/-job[0-9].*$//')
	INPUTFILE2PREFIX=$(echo ${INPUTFILE2##*/} | sed 's/\.fastq.gz$//' | sed 's/\.fq.gz$//' | sed 's/-job[0-9].*$//')

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

	#######################################
	##ACTIONS:

	echo "-----> Processing files: ${INPUTFILE1}, ${INPUTFILE2}"

	cat ${INPUTFILE1} ${INPUTFILE2} > ${OUTPUTFILE}

	##Check for errors.
	if [[ ! -s "${OUTPUTFILE}" ]] ; then
		echo -e "ERROR:\t${INPUTFILE1}\t${INPUTFILE2}"
		echo -e "${INPUTFILE1}\t${INPUTFILE2}" >> ${ERRORFILE}
		echo -e "${INPUTFILE1}\t${INPUTFILE2}\t${OUTPUTFILE}\tERROR" >> ${OUTPUTLIST}
	else
		echo -e "${INPUTFILE1}\t${INPUTFILE2}\t${OUTPUTFILE}\tCLEAR" >> ${OUTPUTLIST}
	fi

done

############################################################################
##SAVE CONTROL FILES:

if [[ ! -s "${ERRORFILE}" ]]; then
	
	##README LOG ENTRY:
printf '%b\n' "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input list: ${INPUTLIST}
Output file: ${OUTPUTLIST}
$(cat ${OUTPUTLIST} | grep -v "^INPUTFILE1" | cut -f3 | sed "s/^/Output file: /g")
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