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
##Align and index paired-end FastQ files to a reference genome using STAR.

##Input $1: Output location.
##Input $2: FASTQ file (read 1, R1) (.fastq.gz or .fq.gz).
##Input $3: FASTQ file (read 2, R2) (.fastq.gz or .fq.gz).
##Input $4: Directory with reference genome STAR-index files.
##Output: Directory containing BAM files (.bam), index files and report files.

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
##    -J fastqPE-BAM-STARalign-${<INPUTFILE>##*/} \
##    fastqPE-BAM-STARalign.sh <OUTPUT LOCATION> <INPUT READ 1 FASTQ FILE> <INPUT READ 2 FASTQ FILE> <REFERENCE GENOME STAR INDEX DIRECTORY>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "fastqPE-BAM-STARalign.sh") 

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
module load star/2.7.11a
module load samtools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE1=$(readlink -f $2)
INPUTFILE2=$(readlink -f $3)
REFGENOMESTARINDEX=$(readlink -f $4)

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

OUTPUTFILEXPREFIX=$(echo "${CONSENSUSFILEPREFIX}.STARalign-job${JOBID}")
OUTPUTDIR=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEXPREFIX}")
OUTPUTFILEPREFIX=$(echo "${OUTPUTDIR}/${OUTPUTFILEXPREFIX}.")

############################################################################
##ACTIONS:

##Create output directory.
mkdir -p ${OUTPUTDIR}

##Align reads.
STAR \
	--readFilesIn ${INPUTFILE1} ${INPUTFILE2} \
	--genomeDir ${REFGENOMESTARINDEX} \
	--outFileNamePrefix ${OUTPUTFILEPREFIX} \
	--readFilesCommand zcat \
	--outFilterType BySJout \
	--outFilterMismatchNoverReadLmax 0.04 \
	--outFilterMultimapNmax 20 \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes Standard \
	--outSAMmapqUnique 255 \
	--twopassMode Basic \
	--chimSegmentMin 15 \
	--chimScoreMin 15 \
	--chimMultimapNmax 5 \
	--chimScoreSeparation 10 \
	--chimJunctionOverhangMin 20 \
	--quantMode TranscriptomeSAM GeneCounts \
	--quantTranscriptomeBan IndelSoftclipSingleend \
	--runThreadN 20

##Index output BAM file.
OUTPUTBAMFILE=$(echo "${OUTPUTFILEPREFIX}Aligned.sortedByCoord.out.bam")

##Sort BAM file.
samtools sort -o ${OUTPUTBAMFILE} ${OUTPUTBAMFILE}

##Index BAM file.
samtools index -b ${OUTPUTBAMFILE}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTBAMFILE}" ]]; then
	
	##README LOG ENTRY:
printf '%b\n' "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE1}
Input file: ${INPUTFILE2}
Input directory: ${REFGENOMESTARINDEX}
Output directory: ${OUTPUTDIR}
" >> $(echo "${OUTPUTLOCATION}/README.txt") 

	##COPY OF SCRIPT:
	cat $0 > $(echo "${OUTPUTLOCATION}/job${JOBID}.sh")
	cat $0 > $(echo "${OUTPUTDIR}/job${JOBID}.sh")
	##COPY OF SLURM FILE:
	cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTLOCATION}/job${JOBID}.out") 
	cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTDIR}/job${JOBID}.out") 

else

	##COPY OF SCRIPT:
	cat $0 > $(echo "${OUTPUTLOCATION}/job${JOBID}.sh")
	cat $0 > $(echo "${OUTPUTDIR}/job${JOBID}.sh")
	##COPY OF SLURM FILE:
	cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTLOCATION}/job${JOBID}.err") 
	cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTDIR}/job${JOBID}.err")

fi

exit 0