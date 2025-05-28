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
##Align and index paired-end FastQ files to a reference genome using BWA.

##Input $1: Output location.
##Input $2: List (.txt) with three tab-separated columns with sample names and FASTQ files (R1\tR2). 
##Input $3: Path to reference genome BWA index file prefix.
##Output: BAM files (.bam) and index files.

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
##    -J list-fastqPE-BAM-BWAalign-${<INPUTFILE>##*/} \
##    list-fastqPE-BAM-BWAalign.sh <OUTPUT LOCATION> <INPUT FILE> <REFERENCE GENOME INDEX FILE PREFIX>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "list-fastqPE-BAM-BWAalign.sh") 

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
module load bcftools
module load bwa

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTLIST=$(readlink -f $2)
GENOMEINDEXFILEPREFIX=$(readlink -f $3)

INPUTLISTLOCATION=${INPUTLIST%/*}
INPUTLISTNAME=${INPUTLIST##*/}

############################################################################
##OUTPUT:

OUTPUTLIST=$(echo "${OUTPUTLOCATION}/files.BWAalign-job${JOBID}.txt")
ERRORFILE=$(echo "${OUTPUTLOCATION}/error.BWAalign-job${JOBID}.txt")

############################################################################
##ACTIONS:

echo -e "SAMPLE\tINPUTFILE1\tINPUTFILE2\tOUTPUTFILE1\tOUTPUTFILE2\tERROR?" > ${OUTPUTLIST}
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

	OUTPUTFILEPREFIX=$(echo "${CONSENSUSFILEPREFIX}.BWAalign-job${JOBID}")
	OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.bam")
	OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.bam.bai")

	TMPFILEPREFIX=$(echo "${OUTPUTLOCATION}/tmp.${OUTPUTFILEPREFIX}")

	#######################################
	##ACTIONS:

	echo "-----> Processing sample: ${SAMPLENAME}"

	##Align reads and sort output BAM file.
	bwa mem -M -t 10 ${GENOMEINDEXFILEPREFIX} ${INPUTFILE1} ${INPUTFILE2} | \
		samtools view -S --bam --with-header - | \
		samtools sort -T ${TMPFILEPREFIX} -o ${OUTPUTFILE1} -
	
	##Index BAM file.
	samtools index --threads 10 -b ${OUTPUTFILE1}

	##Check for errors.
	if [[ ! -s "${OUTPUTFILE1}" ]] || [[ ! -s "${OUTPUTFILE2}" ]] ; then
		echo -e "ERROR:\t${SAMPLENAME}\t${INPUTFILE1}\t${INPUTFILE2}"
		echo -e "${SAMPLENAME}\t${INPUTFILE1}\t${INPUTFILE2}" >> ${ERRORFILE}
		echo -e "${SAMPLENAME}\t${INPUTFILE1}\t${INPUTFILE2}\t${OUTPUTFILE1}\t${OUTPUTFILE2}\tERROR" >> ${OUTPUTLIST}
	else
		echo -e "${SAMPLENAME}\t${INPUTFILE1}\t${INPUTFILE2}\t${OUTPUTFILE1}\t${OUTPUTFILE2}\tCLEAR" >> ${OUTPUTLIST}
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
Input reference genome BWA index file prefix: ${GENOMEINDEXFILEPREFIX}
Output file: ${OUTPUTLIST}
$(cat ${OUTPUTLIST} | grep -v "^SAMPLE" | cut -f4,5 | sed "s/\t/\n/g" | sed "s/^/Output file: /g")
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