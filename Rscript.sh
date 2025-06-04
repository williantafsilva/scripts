#!/bin/bash -l
##Add -l to the shebang to inherit bash profile variables and configuration.
##Required environment variables: ${PATHTOMYSCRIPTS}, ${PATHTOMYSUBMITTEDSCRIPTS}, ${PATHTOMYSLURM}, ${PATHTOPROJTMP}, ${PROJECT_ID}, ${MYSLURMFILE}, ${MYEMAIL}.
############################################################################
################################## SCRIPT ##################################
############################################################################
##Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
##SCRIPT DESCRIPTION:

##Description:
##Run R script (.R).

##Input $1: R script (.R).
##Input $2: Output location.
##Input $3-: R script input arguments.
##Output: R script output.

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
##    -J Rscript-${<INPUTFILE>##*/} \
##    Rscript.sh <R SCRIPT FILE NAME> <OUTPUT LOCATION> <INPUT ARGUMENTS>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "Rscript.sh") 

############################################################################
##JOB ID:

RUNDATE=$(date +"%Y%m%d%H%M%S")
if [[ -z "${SLURM_JOB_ID}" ]] ; then JOBID=${RUNDATE} ; else JOBID=${SLURM_JOB_ID} ; fi 

############################################################################
##INPUT:

RSCRIPTNAME=$1
OUTPUTLOCATION=$(readlink -f $2)
ARGS=${@:2}

RSCRIPT=$(echo "${PATHTOMYSCRIPTS}/${RSCRIPTNAME}") 

############################################################################
##SUBMITTED SCRIPT COPY:

cat $0 > $(echo "${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}") 
cat ${RSCRIPT} > $(echo "${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${RSCRIPTNAME}")

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

module load PDC/24.11
module load R/4.4.2-cpeGNU-24.11

############################################################################
##ACTIONS:

##Set working directory.
cd $(readlink -f .)

##Run R script
IFS=$' '
Rscript --vanilla ${RSCRIPT} ${JOBID} ${ARGS}
IFS=${ORIGINALIFS}

############################################################################
##SAVE CONTROL FILES:

##COPY OF SCRIPTS:
cat $0 > $(echo "${OUTPUTLOCATION}/job${JOBID}.sh")
cat ${RSCRIPT} > $(echo "${OUTPUTLOCATION}/job${JOBID}.R")
##COPY OF SLURM FILE:
cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTLOCATION}/job${JOBID}.out")

exit 0