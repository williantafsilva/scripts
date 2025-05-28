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
##Generate SweepFinder2 output files.

##Input $1: Output location.
##Input $2: SweepFinder2 input file (.SF2input).
##Input $3: Site frequency spectrum file (.SF2sfs).
##Output: SweepFinder2 output file (.SF2output).

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
##    -J SF2input-SF2output-${<INPUTFILE>##*/} \
##    SF2input-SF2output.sh <OUTPUT LOCATION> <SF2INPUT FILE> <SF2SFS FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "SF2input-SF2output.sh") 

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
module load SweepFinder2

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)
INPUTSFSFILE=$(readlink -f $3) ##SFS file.

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.SF2input.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILENAME=$(echo "${OUTPUTFILEPREFIX}.SF2-job${JOBID}.SF2output") 
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}")

############################################################################
##ACTIONS:

##Calculate parameters.
SEQINTERVAL=10000 ##Genetic interval (bp).
POSMAX=$(tail -n+2 ${INPUTFILE} | awk 'NR==1{max = $1 + 0; next} {if ($1 > max) max = $1;} END {print max}') ##Max position.
POSMIN=$(tail -n+2 ${INPUTFILE} | awk 'NR==1{min = $1 + 0; next} {if ($1 < min) min = $1;} END {print min}') ##Min position.
NSITES=$(echo "$(((${POSMAX}-${POSMIN})/${SEQINTERVAL}+1))") ##Number of sites (grid size).

cd ${OUTPUTLOCATION}
if [[ -z "${INPUTSFSFILE}" ]] ; then
    ##Scan for selective sweeps without pre-computed empirical spectrum (.SF2sfs).
    SweepFinder2 -s ${NSITES} ${INPUTFILE} ${OUTPUTFILE}
else
    ##Scan for selective sweeps with pre-computed empirical spectrum (.SF2sfs).
    SweepFinder2 -l ${NSITES} ${INPUTFILE} ${INPUTSFSFILE} ${OUTPUTFILE}
fi

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE}" ]]; then
    
    ##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Input SFS file: ${INPUTSFSFILE}
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