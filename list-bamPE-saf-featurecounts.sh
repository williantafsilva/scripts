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
##Input $2: List (.txt) with paths to sorted and indexed paired-end BAM file (.bam).
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
##    -J list-bamPE-saf-featurecounts-${<INPUTFILE>##*/} \
##    list-bamPE-saf-featurecounts.sh <OUTPUT LOCATION> <INPUT LIST> <INPUT SAF FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "list-bamPE-saf-featurecounts.sh") 

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
INPUTLIST=$(readlink -f $2)
INPUTSAFFILE=$(readlink -f $3)

############################################################################
##OUTPUT:

OUTPUTLIST=$(echo "${OUTPUTLOCATION}/files.featurecounts-job${JOBID}.txt")
ERRORFILE=$(echo "${OUTPUTLOCATION}/error.featurecounts-job${JOBID}.txt")

############################################################################
##ACTIONS:

echo -e "INPUTFILE\tOUTPUTFILE\tERROR?" > ${OUTPUTLIST}
cat ${INPUTLIST} | while read INPUTFILE ; do
	
	#######################################
	##INPUT:

	INPUTFILEPREFIX=$(echo ${INPUTFILE##*/} | sed 's/\.bam$//' | sed 's/-job[0-9].*$//')

	#######################################
	##OUTPUT:

	OUTPUTFILENAME=$(echo "${INPUTFILEPREFIX}.featurecounts-job${JOBID}.txt")
	OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}")

	#######################################
	##ACTIONS:
	
	echo "-----> Procesing file: ${INPUTFILE}"

	##Count features.
	featureCounts -T 10 -F SAF -p --countReadPairs -a ${INPUTSAFFILE} -o ${OUTPUTFILE} ${INPUTFILE}

	##Check for errors.
	if [[ ! -s "${OUTPUTFILE}" ]] ; then
		echo -e "ERROR:\t${INPUTFILE}"
		echo -e "${INPUTFILE}" >> ${ERRORFILE}
		echo -e "${INPUTFILE}\t${OUTPUTFILE}\tERROR" >> ${OUTPUTLIST}
	else
		echo -e "${INPUTFILE}\t${OUTPUTFILE}\tCLEAR" >> ${OUTPUTLIST}
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
Input SAF file: ${INPUTSAFFILE}
Output file: ${OUTPUTLIST}
$(cat ${OUTPUTLIST} | grep -v "^INPUTFILE" | cut -f2 | tr '\t' '\n' | sed "s/^/Output file: /g")
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