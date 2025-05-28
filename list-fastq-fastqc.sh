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
##Run FastQC recursively on all FASTQ files within the input directory.

##Input $1: Output location.
##Input $2: List (.txt) with paths to FASTQ files. 
##Output: HTML files (.html) and ZIP files (.zip).

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
##    -J list-fastq-fastqc-${<INPUTFILE>##*/} \
##    list-fastq-fastqc.sh <OUTPUT LOCATION> <INPUT FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "list-fastq-fastqc.sh") 

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
module load FastQC/0.11.9
module load MultiQC/1.22.2

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTLIST=$(readlink -f $2)

INPUTLISTLOCATION=${INPUTLIST%/*}
INPUTLISTNAME=${INPUTLIST##*/}

############################################################################
##OUTPUT:

OUTPUTLIST=$(echo "${OUTPUTLOCATION}/files.fastqc-job${JOBID}.txt")
ERRORFILE=$(echo "${OUTPUTLOCATION}/error.fastqc-job${JOBID}.txt")

############################################################################
##ACTIONS:

echo -e "INPUTFILE\tOUTPUTFILE1\tOUTPUTFILE2\tERROR?" > ${OUTPUTLIST}
cat ${INPUTLIST} | while read INPUTFILE ; do
	
	#######################################
	##INPUT:

	INPUTFILEPREFIX=$(echo ${INPUTFILE##*/} | sed 's/\.fastq\.gz$//' | sed 's/\.fq\.gz$//')

	#######################################
	##OUTPUT:

	OUTPUTFILEPREFIX=$(echo ${INPUTFILEPREFIX} | sed 's/-job[0-9].*$//')
	OUTPUTFILE1NAME=$(echo "${OUTPUTFILEPREFIX}.fastqc-job${JOBID}_fastqc.html") 
	OUTPUTFILE2NAME=$(echo "${OUTPUTFILEPREFIX}.fastqc-job${JOBID}_fastqc.zip") 
	OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE1NAME}")
	OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE2NAME}")

	#######################################
	##ACTIONS:
	
	echo "-----> Procesing file: ${INPUTFILE}"

	fastqc -o ${OUTPUTLOCATION} ${INPUTFILE}
	mv $(echo "${OUTPUTLOCATION}/${INPUTFILEPREFIX}_fastqc.html") ${OUTPUTFILE1}
	mv $(echo "${OUTPUTLOCATION}/${INPUTFILEPREFIX}_fastqc.zip") ${OUTPUTFILE2}

	##Check for errors.
	if [[ ! -s "${OUTPUTFILE1}" ]] || [[ ! -s "${OUTPUTFILE2}" ]] ; then
		echo -e "ERROR:\t${INPUTFILE}"
		echo -e "${INPUTFILE}" >> ${ERRORFILE}
		echo -e "${INPUTFILE}\t${OUTPUTFILE1}\t${OUTPUTFILE2}\tERROR" >> ${OUTPUTLIST}
	else
		echo -e "${INPUTFILE}\t${OUTPUTFILE1}\t${OUTPUTFILE2}\tCLEAR" >> ${OUTPUTLIST}
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
$(cat ${OUTPUTLIST} | grep -v "^INPUTFILE" | cut -f2,3 | tr '\t' '\n' | sed "s/^/Output file: /g")
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