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
##Calculate coverage based on input annotations file (e.g., BED for genes, transcripts or exons).

##Input $1: Output location.
##Input $2: BED alignment file (.bed) or sorted and indexed BAM file aligned to coordinates (.bam).
##Input $3: Reference annotations file (.bed, .gtf or .gff).
##Output: Coverage data (.txt).

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
##    -J bam-bed-coverage-${<INPUTFILE>##*/} \
##    bam-bed-coverage.sh <OUTPUT LOCATION> <INPUT ALIGNMENT FILE> <INPUT BED/GTF/GFF ANNOTATIONS FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "bam-bed-coverage.sh") 

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
module load bamtools
module load samtools
module load bedtools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTALIGNFILE=$(readlink -f $2)
INPUTANNOTFILE=$(readlink -f $3)

INPUTALIGNFILELOCATION=${INPUTALIGNFILE%/*}
INPUTALIGNFILENAME=${INPUTALIGNFILE##*/}

############################################################################
##OUTPUT:

INPUTFILEPREFIX=$(echo ${INPUTALIGNFILENAME} | sed 's/\.bam$//' | sed 's/\.bed$//' | sed 's/-job[0-9].*$//')
OUTPUTFILENAME=$(echo "${INPUTFILEPREFIX}.bamcoverage-job${JOBID}.txt")
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}") 

############################################################################
##ACTIONS:

##Calculate coverage.
bedtools coverage -a ${INPUTANNOTFILE} -b ${INPUTALIGNFILE} > ${OUTPUTFILE}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input BAM/BED file: ${INPUTALIGNFILE}
Input annotations file: ${INPUTANNOTFILE}
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