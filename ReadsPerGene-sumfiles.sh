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
##Input $2: ReadsPerGene.out.tab file (output from STAR).
##Input $3: ReadsPerGene.out.tab file (output from STAR).
##Output: ReadsPerGene.out.tab file with the sum of the corresponding columns in the input files.

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
##    -J ReadsPerGene-sumfiles-${<INPUTFILE>##*/} \
##    ReadsPerGene-sumfiles.sh <OUTPUT LOCATION> <INPUT FILE 1> <INPUT FILE 2>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "ReadsPerGene-sumfiles.sh") 

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

OUTPUTFILEPREFIX=$(echo ${CONSENSUSFILENAME} | sed 's/\.ReadsPerGene.out.tab$//' | sed 's/-job.*$//')
OUTPUTFILENAME=$(echo "${OUTPUTFILEPREFIX}.ReadsPerGeneSum-job${JOBID}.ReadsPerGene.out.tab")
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}")

############################################################################
##ACTIONS:

paste "${INPUTFILE1}" "${INPUTFILE2}" | while IFS=$'\t' read -r F1_C1 F1_C2 F1_C3 F1_C4 F2_C1 F2_C2 F2_C3 F2_C4 ; do
	if [[ "${F1_C1}" == "${F2_C1}" ]]; then
    	##Compute the sum of the columns.
    	SUM_C2=$((F1_C2 + F2_C2))
    	SUM_C3=$((F1_C3 + F2_C3))
    	SUM_C4=$((F1_C4 + F2_C4))
    
    	##Output the row with the sums
    	echo -e "${F1_C1}\t${SUM_C2}\t${SUM_C3}\t${SUM_C4}" >> "${OUTPUTFILE}"
    else 
    	echo "ERROR: Row names (${F1_C1}; ${F2_C1}) do not match."
    	echo -e "ERROR: Row names do not match.\t?\t?\t?" >> "${OUTPUTFILE}"
    fi
done

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE1}
Input file: ${INPUTFILE2}
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