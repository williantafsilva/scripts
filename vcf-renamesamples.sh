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
##Change sample names in a VCF file (.vcf or .vcf.gz) and produce index file (.tbi).

##Input $1: Output location.
##Input $2: VCF file (.vcf or .vcf.gz).
##Input $3: Sample name conversion file (OLDNAME NEWNAME\n).
##Output: VCF file (.vcf.gz) and index file (.tbi).

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
##    -J vcf-renamesamples-${<INPUTFILE>##*/} \
##    vcf-renamesamples.sh <OUTPUT LOCATION> <INPUT VCF FILE> <SAMPLE NAME CONVERTION FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-renamesamples.sh") 

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
module load bcftools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)
SAMPLENAMECONVERTIONFILE=$(readlink -f $3) ##Sample name conversion file.

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/\.bcf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILE1NAME=$(echo "${OUTPUTFILEPREFIX}.renamesamples-job${JOBID}.vcf.gz") 
OUTPUTFILE2NAME=$(echo "${OUTPUTFILEPREFIX}.renamesamples-job${JOBID}.vcf.gz.tbi") 
OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE1NAME}") 
OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE2NAME}") 

############################################################################
##ACTIONS:

bcftools reheader --samples ${SAMPLENAMECONVERTIONFILE} --output ${OUTPUTFILE1} ${INPUTFILE}
bcftools index --tbi --output ${OUTPUTFILE2} ${OUTPUTFILE1}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE1}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Input sample name convertion file: ${SAMPLENAMECONVERTIONFILE}
Output file: ${OUTPUTFILE1}
Output file: ${OUTPUTFILE2}
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