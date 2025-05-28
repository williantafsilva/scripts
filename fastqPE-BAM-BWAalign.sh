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
##Input $2: FASTQ file (read 1, R1) (.fastq.gz or .fq.gz).
##Input $3: FASTQ file (read 2, R2) (.fastq.gz or .fq.gz).
##Input $4: Path to file prefix of BWA index of reference genome FASTA file.
##Output: BAM file (.bam) and index file.

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
##    -J fastqPE-BAM-BWAalign-${<INPUTFILE>##*/} \
##    fastqPE-BAM-BWAalign.sh <OUTPUT LOCATION> <INPUT READ 1 FASTQ FILE> <INPUT READ 2 FASTQ FILE> <REFERENCE GENOME INDEX FILE PREFIX>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "fastqPE-BAM-BWAalign.sh") 

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
INPUTFILE1=$(readlink -f $2)
INPUTFILE2=$(readlink -f $3)
REFGENOMEFILE=$(readlink -f $4)

############################################################################
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

TMPFILEPREFIX=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}")

############################################################################
##ACTIONS:

##Align reads and sort output BAM file.
bwa mem -M -t 10 ${REFGENOMEFILE} ${INPUTFILE1} ${INPUTFILE2} | \
	samtools view -S --bam --with-header - | \
	samtools sort -T ${TMPFILEPREFIX} -o ${OUTPUTFILE1} -
	
##Index BAM file.
samtools index --threads 10 -b ${OUTPUTFILE1}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE1}" ]]; then
	
	##README LOG ENTRY:
printf '%b\n' "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE1}
Input file: ${INPUTFILE2}
Input file: ${REFGENOMEFILE}
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