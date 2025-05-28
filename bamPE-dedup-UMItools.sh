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
##Index BCF file and count variants.

##Input $1: Output location.
##Input $2: Indexed BAM file (*Aligned.sortedByCoord.out.bam).
##Output: Deduplicated BAM file (.bam), index file (.bai), log file (.log) and stats files (.tsv).

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
##    -J bamPE-dedup-UMItools-${<INPUTFILE>##*/} \
##    bamPE-dedup-UMItools.sh <OUTPUT LOCATION> <INPUT BAM FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "bamPE-dedup-UMItools.sh") 

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
module load umi_tools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.bam.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILEPREFIX=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.UMIdedup-job${JOBID}")
OUTPUTFILE1=$(echo "${OUTPUTFILEPREFIX}.bam") 
OUTPUTFILE2=$(echo "${OUTPUTFILEPREFIX}.log") 

############################################################################
##ACTIONS:

umi_tools dedup \
	-I ${INPUTFILE} \
	-S ${OUTPUTFILE1} \
	-L ${OUTPUTFILE2} \
	--paired \
	--buffer-whole-contig \
	--multimapping-detection-method=NH \
	--output-stats=${OUTPUTFILEPREFIX} \
	--umi-separator=':'

##Sort BAM file.
samtools sort -o ${OUTPUTFILE1} ${OUTPUTFILE1}

##Index BAM file.
samtools index -b ${OUTPUTFILE1}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE1}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Output file: ${OUTPUTFILE1}
Output file: ${OUTPUTFILE1}.bai
Output file: ${OUTPUTFILE2}
Output file: ${OUTPUTFILEPREFIX}*.tsv
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