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
##Retrieve sequence from indexed reference FASTA file.

##Input $1: Output location.
##Input $2: Sequence range (chr:position-position).
##Input $3: Reference FASTA file (.fa).
##Input $4: Target FASTA file (.fa).
##Input $5: File tag (to be included in the output file name).
##Output: FASTA file (.fa).

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
##    -J blastn-seqrange-refFASTA-targetFASTA-${<INPUTFILE>##*/} \
##    blastn-seqrange-refFASTA-targetFASTA.sh <OUTPUT LOCATION> <SEQUENCE RANGE> <INPUT REFERENCE FASTA FILE> <INPUT TARGET FASTA FILE> <FILE TAG>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "blastn-seqrange-refFASTA-targetFASTA.sh") 

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
module load blast/2.15.0+

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
SEQRANGE=$2 #Sequence range.
INPUTREFFASTAFILE=$(readlink -f $3)
INPUTTARGETFASTAFILE=$(readlink -f $4)
FILETAG=$5 #File tag.

INPUTREFFASTAFILELOCATION=${INPUTREFFASTAFILE%/*}
INPUTREFFASTAFILENAME=${INPUTREFFASTAFILE##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTREFFASTAFILENAME} | sed 's/\.fa.*$//' | sed 's/\.fasta.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILE1NAME=$(echo "${OUTPUTFILEPREFIX}.blastn${FILETAG}-job${JOBID}.queryseq.fa")
OUTPUTFILE2NAME=$(echo "${OUTPUTFILEPREFIX}.blastn${FILETAG}-job${JOBID}.tophits.txt")
OUTPUTFILE3NAME=$(echo "${OUTPUTFILEPREFIX}.blastn${FILETAG}-job${JOBID}.alignments.txt")
OUTPUTFILE4NAME=$(echo "${OUTPUTFILEPREFIX}.blastn${FILETAG}-job${JOBID}.tophitseq.fa")
OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE1NAME}") 
OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE2NAME}") 
OUTPUTFILE3=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE3NAME}") 
OUTPUTFILE4=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE4NAME}") 

#Temporary files.
TMPFILE1=$(echo "${OUTPUTLOCATION}/tmp-blasttargetdb-job${JOBID}") 

############################################################################
##ACTIONS:

#Index reference FASTA file.
samtools faidx "${INPUTREFFASTAFILE}"

#Extract the query sequence.
samtools faidx "${INPUTREFFASTAFILE}" "${SEQRANGE}" > "${OUTPUTFILE1}"

#Create BLAST database.
makeblastdb \
    -in "${INPUTTARGETFASTAFILE}" \
    -dbtype nucl \
    -out "${TMPFILE1}" \
    >/dev/null

#Run BLAST and report top 10 hits in tabular format.
#Coordinates are retained in their original orientation.
#-max_target_seqs limits the number of different subject/database sequences (target hits) reported or kept.
#-max_hsps limits the number of High-scoring Segment Pairs (HSPs / local alignments) allowed for each individual query-subject pair.
blastn \
    -query "${OUTPUTFILE1}" \
    -db "${TMPFILE1}" \
    -task blastn \
    -outfmt "6 qseqid qlen sseqid slen qstart qend sstart send nident pident length mismatch gaps qcovs evalue bitscore" \
    -max_target_seqs 10 \
    -max_hsps 1 \
    -out "${OUTPUTFILE2}"
#-evalue 1e-05 \
#-perc_identity 60 \
#-num_threads 10 \

#Append column names.
sed -i '1i qseqid\tqlen\tsseqid\tslen\tqstart\tqend\tsstart\tsend\tnident\tpident\tlength\tmismatch\tgaps\tqcovs\tevalue\tbitscore' "${OUTPUTFILE2}"

#Print BLASTN results.
echo -e "#BLASTN output:\n"
cat "${OUTPUTFILE2}"

#Run BLAST and report alignments of top 10 hits in pairwise format.
blastn \
    -query "${OUTPUTFILE1}" \
    -db "${TMPFILE1}" \
    -task blastn \
    -outfmt 0 \
    -max_target_seqs 10 \
    -max_hsps 1 \
    -out "${OUTPUTFILE3}"

#Select the first hit reported by BLAST.
FIRSTHIT=$(cat "${OUTPUTFILE2}" | tail -n+2 | head -n1)
if [[ ! -z "${FIRSTHIT}" ]]; then
    #Split the first hit into different variables.
	TARGET_CHR=$(cat "${OUTPUTFILE2}" | tail -n+2 | head -n1 | cut -f3)
	TARGET_START=$(cat "${OUTPUTFILE2}" | tail -n+2 | head -n1 | cut -f7)
	TARGET_END=$(cat "${OUTPUTFILE2}" | tail -n+2 | head -n1 | cut -f8)
	#samtools faidx expects ascending coordinates, regardless of strand.
	if (( TARGET_START <= TARGET_END )); then
	    TARGET_REGION="${TARGET_CHR}:${TARGET_START}-${TARGET_END}"
	    STRAND="+"
	else
	    TARGET_REGION="${TARGET_CHR}:${TARGET_END}-${TARGET_START}"
	    STRAND="-"
	fi
	#Extract the corresponding sequence from the target genome.
	#samtools faidx returns the sequence in reference-forward orientation.
	samtools faidx --output "${OUTPUTFILE4}" "${INPUTTARGETFASTAFILE}" "${TARGET_REGION}" 
else
	echo "No BLAST hits found for ${SEQRANGE}."
	touch "${OUTPUTFILE4}"
fi

#Delete temporary files.
rm -f "${TMPFILE1}"*

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE2}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Sequence range: ${SEQRANGE}
Input reference FASTA file: ${INPUTREFFASTAFILE}
Input target FASTA file: ${INPUTTARGETFASTAFILE}
File tag: ${FILETAG}
Output file: ${OUTPUTFILE1}
Output file: ${OUTPUTFILE2}
Output file: ${OUTPUTFILE3}
Output file: ${OUTPUTFILE4}
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