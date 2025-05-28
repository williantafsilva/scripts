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
##Pre-process paired-end FastQ files using Fastp.

##Input $1: Output location.
##Input $2: FASTQ file (read 1, R1) (.fastq.gz or .fq.gz).
##Input $3: FASTQ file (read 2, R2) (.fastq.gz or .fq.gz).
##Output: FASTQ files (.fastq.gz) and report files (.json and .html).

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
##    -J fastqPE-fastp-${<INPUTFILE>##*/} \
##    fastqPE-fastp.sh <OUTPUT LOCATION> <INPUT READ 1 FASTQ FILE> <INPUT READ 2 FASTQ FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "fastqPE-fastp.sh") 

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
module load fastp/0.23.4

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

##Create consensus file name.
CONSENSUSFILENAME=""
while read C ; do
    C1=$(echo ${C} | cut -f1)
    C2=$(echo ${C} | cut -f2)
    if [[ "${C1}" == "${C2}" ]] ; then
        CONSENSUSFILENAME+="${C1}"
    else
        CONSENSUSFILENAME+="X"
    fi
done <<< "$(paste -d'\t' <(echo ${INPUTFILE1NAME} | grep -o .) <(echo ${INPUTFILE2NAME} | grep -o .))"

OUTPUTFILE1PREFIX=$(echo ${INPUTFILE1NAME} | sed 's/\.fastq\..*$//' | sed 's/\.fq\..*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILE2PREFIX=$(echo ${INPUTFILE2NAME} | sed 's/\.fastq\..*$//' | sed 's/\.fq\..*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILEXPREFIX=$(echo ${CONSENSUSFILENAME} | sed 's/\.fastq\..*$//' | sed 's/\.fq\..*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILE1NAME=$(echo "${OUTPUTFILE1PREFIX}.fastp-job${JOBID}.fastq.gz") 
OUTPUTFILE2NAME=$(echo "${OUTPUTFILE2PREFIX}.fastp-job${JOBID}.fastq.gz") 
OUTPUTFILE3NAME=$(echo "${OUTPUTFILEXPREFIX}.fastp_report-job${JOBID}.json") 
OUTPUTFILE4NAME=$(echo "${OUTPUTFILEXPREFIX}.fastp_report-job${JOBID}.html") 
OUTPUTFILE5NAME=$(echo "${OUTPUTFILE1PREFIX}.fastp_failQC_unpaired1-job${JOBID}.fastq.gz") 
OUTPUTFILE6NAME=$(echo "${OUTPUTFILE2PREFIX}.fastp_failQC_unpaired2-job${JOBID}.fastq.gz") 
OUTPUTFILE7NAME=$(echo "${OUTPUTFILEXPREFIX}.fastp_failQC_filters-job${JOBID}.fastq.gz") 
OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE1NAME}") 
OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE2NAME}") 
OUTPUTFILE3=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE3NAME}") 
OUTPUTFILE4=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE4NAME}") 

##Create output directory for reads that fail quality control.
FAILQCDIR=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEXPREFIX}.fastp_failQCreads-job${JOBID}") 
mkdir -p ${FAILQCDIR}

OUTPUTFILE5=$(echo "${FAILQCDIR}/${OUTPUTFILE5NAME}") 
OUTPUTFILE6=$(echo "${FAILQCDIR}/${OUTPUTFILE6NAME}") 
OUTPUTFILE7=$(echo "${FAILQCDIR}/${OUTPUTFILE7NAME}") 

############################################################################
##ACTIONS:

fastp \
	--in1 ${INPUTFILE1} \
	--in2 ${INPUTFILE2} \
	--out1 ${OUTPUTFILE1} \
	--out2 ${OUTPUTFILE2} \
	--json ${OUTPUTFILE3} \
	--html ${OUTPUTFILE4} \
	--unpaired1 ${OUTPUTFILE5} \
	--unpaired2 ${OUTPUTFILE6} \
	--failed_out ${OUTPUTFILE7} \
	--dont_overwrite \
	--trim_poly_x \
	--trim_poly_g \
	--umi \
	--umi_loc read1 \
	--umi_len 12 \
	--dedup \
	--overrepresentation_analysis \
	--overrepresentation_sampling 100 \
	--detect_adapter_for_pe \
	--qualified_quality_phred 20 \
	--thread 10

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE1}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE1}
Input file: ${INPUTFILE2}
Output file: ${OUTPUTFILE1}
Output file: ${OUTPUTFILE2}
Output file: ${OUTPUTFILE3}
Output file: ${OUTPUTFILE4}
Output directory: ${FAILQCDIR}
Output file: ${OUTPUTFILE5}
Output file: ${OUTPUTFILE6}
Output file: ${OUTPUTFILE7}
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