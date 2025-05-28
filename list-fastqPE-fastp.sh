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
##Input $2: List (.txt) with three tab-separated columns with sample names and FASTQ files (R1\tR2). 
##Output: Directories with FASTQ files (.fastq.gz) and report files (.json and .html) for each sample.

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
##    -J list-fastqPE-fastp-${<INPUTFILE>##*/} \
##    list-fastqPE-fastp.sh <OUTPUT LOCATION> <INPUT FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "list-fastqPE-fastp.sh") 

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
INPUTLIST=$(readlink -f $2)

INPUTLISTLOCATION=${INPUTLIST%/*}
INPUTLISTNAME=${INPUTLIST##*/}

############################################################################
##OUTPUT:

OUTPUTLIST=$(echo "${OUTPUTLOCATION}/files.fastp-job${JOBID}.txt")
ERRORFILE=$(echo "${OUTPUTLOCATION}/error.fastp-job${JOBID}.txt")

############################################################################
##ACTIONS:

echo -e "SAMPLE\tINPUTFILE1\tINPUTFILE2\tOUTPUTDIR\tOUTPUTFILE1\tOUTPUTFILE2\tOUTPUTFILE3\tOUTPUTFILE4\tOUTPUTFILE5\tOUTPUTFILE6\tOUTPUTFILE7\tERROR?" > ${OUTPUTLIST}
cat ${INPUTLIST} | while read LINE ; do
	
	#######################################
	##INPUT:

	SAMPLENAME=$(echo ${LINE} | cut -f1)
	INPUTFILE1=$(echo ${LINE} | cut -f2)
	INPUTFILE2=$(echo ${LINE} | cut -f3)

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

	##Create sample output directory.
	OUTPUTDIR=$(echo "${OUTPUTLOCATION}/${CONSENSUSFILEPREFIX}.fastp-job${JOBID}") 
	mkdir -p ${OUTPUTDIR}

	OUTPUTFILE1NAME=$(echo "${INPUTFILE1PREFIX}.fastp-job${JOBID}.fastq.gz") 
	OUTPUTFILE2NAME=$(echo "${INPUTFILE2PREFIX}.fastp-job${JOBID}.fastq.gz") 
	OUTPUTFILE3NAME=$(echo "${CONSENSUSFILEPREFIX}.fastp_report-job${JOBID}.json") 
	OUTPUTFILE4NAME=$(echo "${CONSENSUSFILEPREFIX}.fastp_report-job${JOBID}.html") 
	OUTPUTFILE5NAME=$(echo "${INPUTFILE1PREFIX}.fastp_failQC_unpaired1-job${JOBID}.fastq.gz") 
	OUTPUTFILE6NAME=$(echo "${INPUTFILE2PREFIX}.fastp_failQC_unpaired2-job${JOBID}.fastq.gz") 
	OUTPUTFILE7NAME=$(echo "${CONSENSUSFILEPREFIX}.fastp_failQC_filters-job${JOBID}.fastq.gz") 
	OUTPUTFILE1=$(echo "${OUTPUTDIR}/${OUTPUTFILE1NAME}") 
	OUTPUTFILE2=$(echo "${OUTPUTDIR}/${OUTPUTFILE2NAME}") 
	OUTPUTFILE3=$(echo "${OUTPUTDIR}/${OUTPUTFILE3NAME}") 
	OUTPUTFILE4=$(echo "${OUTPUTDIR}/${OUTPUTFILE4NAME}") 
	OUTPUTFILE5=$(echo "${OUTPUTDIR}/${OUTPUTFILE5NAME}") 
	OUTPUTFILE6=$(echo "${OUTPUTDIR}/${OUTPUTFILE6NAME}") 
	OUTPUTFILE7=$(echo "${OUTPUTDIR}/${OUTPUTFILE7NAME}") 

	#######################################
	##ACTIONS:

	echo "-----> Processing sample: ${SAMPLENAME}"

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

	##Check for errors.
	if [[ ! -s "${OUTPUTFILE1}" ]] || [[ ! -s "${OUTPUTFILE2}" ]] || [[ ! -s "${OUTPUTFILE3}" ]] || [[ ! -s "${OUTPUTFILE4}" ]] || [[ ! -s "${OUTPUTFILE5}" ]] || [[ ! -s "${OUTPUTFILE6}" ]] || [[ ! -s "${OUTPUTFILE7}" ]] ; then
		echo -e "ERROR:\t${SAMPLENAME}\t${INPUTFILE1}\t${INPUTFILE2}"
		echo -e "${SAMPLENAME}\t${INPUTFILE1}\t${INPUTFILE2}" >> ${ERRORFILE}
		echo -e "${SAMPLENAME}\t${INPUTFILE1}\t${INPUTFILE2}\t${OUTPUTDIR}\t${OUTPUTFILE1}\t${OUTPUTFILE2}\t${OUTPUTFILE3}\t${OUTPUTFILE4}\t${OUTPUTFILE5}\t${OUTPUTFILE6}\t${OUTPUTFILE7}\tERROR" >> ${OUTPUTLIST}
	else
		echo -e "${SAMPLENAME}\t${INPUTFILE1}\t${INPUTFILE2}\t${OUTPUTDIR}\t${OUTPUTFILE1}\t${OUTPUTFILE2}\t${OUTPUTFILE3}\t${OUTPUTFILE4}\t${OUTPUTFILE5}\t${OUTPUTFILE6}\t${OUTPUTFILE7}\tCLEAR" >> ${OUTPUTLIST}
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
$(cat ${OUTPUTLIST} | grep -v "^SAMPLE" | cut -f4 | sed "s/^/Output directory: /g")
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

cat ${OUTPUTLIST} | grep -v "^SAMPLE" | grep -vP "\tERROR" | while read L ; do
	
	OUTPUTDIR=$(echo ${L} | cut -f4)

	##README LOG ENTRY:
printf '%b\n' "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input list: ${INPUTLIST}
Output file: ${OUTPUTLIST}
Output directory: $(echo -e ${L} | cut -f4 | sed "s/^/Output directory: /g")
$(echo -e ${L} | cut -f5-11 | tr '\t' '\n' | sed "s/^/Output file: /g")
" >> $(echo "${OUTPUTDIR}/README.txt") 

	##COPY OF SCRIPT:
	cat $0 > $(echo "${OUTPUTDIR}/job${JOBID}.sh")
	##COPY OF SLURM FILE:
	cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTDIR}/job${JOBID}.out") 

done

exit 0