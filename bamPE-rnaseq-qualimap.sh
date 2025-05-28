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
##Merge BAM files (.bam).

##Input $1: Output location.
##Input $2: Paired-end RNAseq BAM file (.bam).
##Input $3: Reference genome annotations file (.gtf).
##Output: Directory with quality control report and data.

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
##    -J bamPE-rnaseq-qualimap-${<INPUTFILE>##*/} \
##    bamPE-rnaseq-qualimap.sh <OUTPUT LOCATION> <INPUT BAM FILE> <INPUT GTF FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "bamPE-rnaseq-qualimap.sh") 

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

module load PDC/23.12 
module load bioinfo-tools
module load QualiMap/2.2.1

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTBAMFILE=$(readlink -f $2)
INPUTGTFFILE=$(readlink -f $3)

INPUTBAMFILELOCATION=${INPUTBAMFILE%/*}
INPUTBAMFILENAME=${INPUTBAMFILE##*/}

############################################################################
##OUTPUT:

INPUTBAMFILEPREFIX=$(echo ${INPUTBAMFILENAME} | sed 's/\.bam.*$//' | sed 's/-job[0-9].*$//')
OUTPUTDIRNAME=$(echo "${INPUTBAMFILEPREFIX}.qualimap-job${JOBID}")
OUTPUTDIR=$(echo "${OUTPUTLOCATION}/${OUTPUTDIRNAME}")

############################################################################
##ACTIONS:

##Create output directory.
mkdir -p ${OUTPUTDIR}

qualimap rnaseq --paired --algorithm uniquely-mapped-reads -bam ${INPUTBAMFILE} -gtf ${INPUTGTFFILE} -outdir ${OUTPUTDIR}

############################################################################
##SAVE CONTROL FILES:

if [[ $(ls ${OUTPUTDIR} | wc -l) -gt 0 ]]; then
    
    ##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTBAMFILE}
Input file: ${INPUTGTFFILE}
Output directory: ${OUTPUTDIR}
" >> $(echo "${OUTPUTLOCATION}/README.txt") 

echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTBAMFILE}
Input file: ${INPUTGTFFILE}
Output directory: ${OUTPUTDIR}
" >> $(echo "${OUTPUTDIR}/README.txt") 

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