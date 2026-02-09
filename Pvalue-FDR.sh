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

#Create a temporary files.
TMP1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.PvalueFDR-tmp1-job${JOBID}.txt") 
TMP2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.PvalueFDR-tmp2-job${JOBID}.txt") 
TMP3=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.PvalueFDR-tmp3-job${JOBID}.txt") 
TMP4=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.PvalueFDR-tmp4-job${JOBID}.txt") 

############################################################################
##ACTIONS:

##Create a temporary file.
#TMP=$(mktemp)
#
##Delete temporary file on exit (normal or error).
#trap 'rm -f "${TMP}"' EXIT
#
#Process data.
#tail -n +2 "${INPUTFILE}" \
#| awk -v PCOL="${PVALUECOLUMNINDEX}" 'BEGIN{OFS="\t"} {print NR, $PCOL, $0}' \
#| sort -g -k2,2 \
#| awk -v NTESTS="${NTESTS}" 'BEGIN{OFS="\t"}
#{
#    rank = NR #BH rank (1 = smallest p-value).
#    idx[rank] = $1 #Map rank -> original row index.
#    line[rank] = substr($0, index($0,$3)) #Store original line (everything from column 3 onward).
#    bh[rank] = $2 * NTESTS / rank #Benjamini–Hochberg formula.
#}
#END {
#    #Enforce monotonicity. Adjusted p-values must be non-decreasing.
#    for (i = NR - 1; i >= 1; i--) {
#        if (bh[i] > bh[i + 1])
#            bh[i] = bh[i + 1]
#    }
#
#    #Cap FDR p-values at 1 and restore original order.
#    for (i = 1; i <= NR; i++) {
#        adj[idx[i]] = (bh[i] > 1 ? 1 : bh[i])
#    }
#
#    #Print lines in original order with appended FDR value.
#    for (i = 1; i <= NR; i++) {
#        printf "%s\t%.10g\n", line[i], adj[i]
#    }
#}' > "${TMP}"
#
##Save header.
#awk 'NR==1 {print $0 "\tPvalueFDR"}' "${INPUTFILE}" > ${OUTPUTFILE}
#
##Save data.
#cat "${TMP}" >> ${OUTPUTFILE}

#Select target column containing raw p-values.
cat ${INPUTFILE} | cut -f ${PVALUECOLUMNINDEX} > ${TMP1} 

#Add original order of input p-values, sort and add column with order of ascending raw p-values.
cat ${TMP1} | awk -v OFS='\t' 'NR==1 { print $0, "OriginalOrder"; next } { print $0, NR-1 }' \
| tail -n+2 \
| sort -g -k1,1 \
| awk -v OFS='\t' '{ print $0, NR }' > ${TMP2}

#Calculate BH-FDR p-values.
cat ${TMP2} | awk -v OFS='\t' -v NTESTS=${NTESTS} '
{
    p[NR] = $1
    rank[NR] = $3
    line[NR] = $0
    q[NR] = p[NR] * NTESTS / rank[NR]
}
END {
    #Enforce monotonicity from bottom to top.
    for (i = NR - 1; i >= 1; i--) {
        if (q[i] > q[i+1])
            q[i] = q[i+1]
    }

    #Print results.
    for (i = 1; i <= NR; i++) {
        if (q[i] > 1) q[i] = 1
        print line[i], q[i]
    }
}
' > ${TMP3} 

#Sort p-values according to original order.
cat ${TMP3} | sort -g -k2,2 | sed '1ip-value\tOriginalOrder\tAscendingOrder\tPvalueFDR' > ${TMP4}

#Save data.
paste ${INPUTFILE} <(cut -f4 ${TMP4}) > ${OUTPUTFILE}

#Delete temporary file.
rm -f ${TMP1} ${TMP2} ${TMP3} ${TMP4}

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