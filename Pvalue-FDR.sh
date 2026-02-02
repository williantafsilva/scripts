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
##Calculate genome-wide statistics and coverage of a BAM file.

##Input $1: Output location.
##Input $2: Tab-separated file with header (.txt).
##Input $3: Index of column containing raw p-values.
##Input $4: Total number of tests.
##Output: Input file with an extra column containg FDR p-values.

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
##    -J Pvalue-FDR-${<INPUTFILE>##*/} \
##    Pvalue-FDR.sh <OUTPUT LOCATION> <INPUT FILE> <P-VALUE COLUMN INDEX> <NUMBER OF TESTS>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "Pvalue-FDR.sh") 

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

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)
PVALUECOLUMNINDEX=$3
NTESTS=$4

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.txt.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILENAME=$(echo "${OUTPUTFILEPREFIX}.PvalueFDR-job${JOBID}.txt") 
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}") 

############################################################################
##ACTIONS:

#Create a temporary file.
TMP=$(mktemp)

#Delete temporary file on exit (normal or error).
trap 'rm -f "${TMP}"' EXIT

#Process data.
tail -n +2 "${INPUTFILE}" \
| awk -v PCOL="${PVALUECOLUMNINDEX}" 'BEGIN{OFS="\t"} {print NR, $PCOL, $0}' \
| sort -k2,2g \
| awk -v NTESTS="${NTESTS}" 'BEGIN{OFS="\t"}
{
    rank = NR #BH rank (1 = smallest p-value).
    idx[rank] = $1 #Map rank -> original row index.
    line[rank] = substr($0, index($0,$3)) #Store original line (everything from column 3 onward).
    bh[rank] = $2 * NTESTS / rank #Benjamini–Hochberg formula.
}
END {
    #Enforce monotonicity. Adjusted p-values must be non-decreasing.
    for (i = NR - 1; i >= 1; i--) {
        if (bh[i] > bh[i + 1])
            bh[i] = bh[i + 1]
    }

    #Cap FDR p-values at 1 and restore original order.
    for (i = 1; i <= NR; i++) {
        adj[idx[i]] = (bh[i] > 1 ? 1 : bh[i])
    }

    #Print lines in original order with appended FDR value.
    for (i = 1; i <= NR; i++) {
        printf "%s\t%.10g\n", line[i], adj[i]
    }
}' > "${TMP}"

#Save header.
awk 'NR==1 {print $0 "\tPvalueFDR"}' "${INPUTFILE}" > ${OUTPUTFILE}

#Save data.
cat "${TMP}" >> ${OUTPUTFILE}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Input p-value column index: ${PVALUECOLUMNINDEX}
Input number of tests: ${NTESTS}
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