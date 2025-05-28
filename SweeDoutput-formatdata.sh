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
##Concatenate and format SweeD output data to be used for plotting (.txt).

##Input $1: Output location.
##Input $2: Directory containing SweeD output files (per chromosome).
##Output: Concatenated and formatted output data to be used for plotting (.txt).

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
##    -J SweeDoutput-formatdata-${<INPUTDIR>##*/} \
##    SweeDoutput-formatdata.sh <OUTPUT LOCATION> <INPUT DIRECTORY> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "SweeDoutput-formatdata.sh") 

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
INPUTDIR=$(readlink -f $2)

INPUTDIRLOCATION=${INPUTDIR%/*}
INPUTDIRNAME=${INPUTDIR##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTDIRNAME} | sed 's/-job[0-9].*$//')
OUTPUTFILENAME=$(echo "${OUTPUTFILEPREFIX}.SweeDdata-job${JOBID}.txt") 
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}") 

TMPFILE=$(echo "${OUTPUTLOCATION}/tmp-job${JOBID}.txt") 

############################################################################
##ACTIONS:

##Concatenate output files and add chromosome column.
find ${INPUTDIR}/SweeD_Report.* -maxdepth 0 | while read F ; do
	CHROM=$(echo "${F##*/}" | cut -d. -f2)
	tail -n+2 ${F} | sed "s/^/${CHROM}\t/" >> ${TMPFILE}
done

##Sort output file and add header.
cat ${TMPFILE} | sort -V -k1,1 -k2,2n > ${OUTPUTFILE}
sed -i '1iChromosome\tPosition\tLikelihood\tAlpha\tStartPos\tEndPos' ${OUTPUTFILE}

##Remove temporary files.
rm -f ${TMPFILE}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input directory: ${INPUTDIR}
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