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
##Convert BCF file into compressed VCF file (.vcf.gz).

##Input $1: Output location.
##Input $2: Ensembl/GENCODE GTF file (.gtf).
##Output: SAF files (.saf) for genes, transcripts, exons (with exon ID) and exons (with gene ID).

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
##    -J gtf-saf-${<INPUTFILE>##*/} \
##    gtf-saf.sh <OUTPUT LOCATION> <INPUT GTF FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "gtf-saf.sh") 

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
module load bcftools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.gtf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILE1NAME=$(echo "${OUTPUTFILEPREFIX}.gtftosafgenes-job${JOBID}.saf") 
OUTPUTFILE2NAME=$(echo "${OUTPUTFILEPREFIX}.gtftosaftranscripts-job${JOBID}.saf") 
OUTPUTFILE3NAME=$(echo "${OUTPUTFILEPREFIX}.gtftosafexons-job${JOBID}.saf") 
OUTPUTFILE4NAME=$(echo "${OUTPUTFILEPREFIX}.gtftosafexonsgeneid-job${JOBID}.saf") 
OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE1NAME}") 
OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE2NAME}") 
OUTPUTFILE3=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE3NAME}") 
OUTPUTFILE4=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE4NAME}") 

############################################################################
##ACTIONS:

##Convert GTF file to simplified annotations format (SAF) for genes.
paste -d' ' <(grep -P "\tgene\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
awk -F $'\t' 'BEGIN { OFS=FS } { print $6,$1,$2,$3,$4 }') \
<(grep -P "\tgene\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | sed 's/.*\(ENSGALG[0-9]*\t\).*/\1/g') | \
awk '{print $6"\t"$2"\t"$3"\t"$4"\t"$5}' | \
sort -V -k2,2 -k3,3n -k4,4n > ${OUTPUTFILE1}
sed -i '1iGeneID\tChr\tStart\tEnd\tStrand' ${OUTPUTFILE1}

##Convert GTF file to simplified annotations format (SAF) for transcripts.
paste -d' ' <(grep -P "\ttranscript\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
awk -F $'\t' 'BEGIN { OFS=FS } { print $10,$1,$2,$3,$4 }') \
<(grep -P "\ttranscript\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | sed 's/.*\(ENSGALT[0-9]*\t\).*/\1/g') | \
awk '{print $6"\t"$2"\t"$3"\t"$4"\t"$5}' | \
sort -V -k2,2 -k3,3n -k4,4n > ${OUTPUTFILE2}
sed -i '1iTranscriptID\tChr\tStart\tEnd\tStrand' ${OUTPUTFILE2}

##Convert GTF file to simplified annotations format (SAF) for exons (with exon ID).
paste -d' ' <(grep -P "\texon\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
awk -F $'\t' 'BEGIN { OFS=FS } { print $6,$1,$2,$3,$4 }') \
<(grep -P "\texon\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | sed 's/.*\(ENSGALE[0-9]*\t\).*/\1/g') | \
awk '{print $6"\t"$2"\t"$3"\t"$4"\t"$5}' | \
sort -V -k2,2 -k3,3n -k4,4n > ${OUTPUTFILE3}
sed -i '1iExonID\tChr\tStart\tEnd\tStrand' ${OUTPUTFILE3}

##Convert GTF file to simplified annotations format (SAF) for exons (with gene ID).
paste -d' ' <(grep -P "\texon\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
awk -F $'\t' 'BEGIN { OFS=FS } { print $6,$1,$2,$3,$4 }') \
<(grep -P "\texon\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | sed 's/.*\(ENSGALG[0-9]*\t\).*/\1/g') | \
awk '{print $6"\t"$2"\t"$3"\t"$4"\t"$5}' | \
sort -V -k2,2 -k3,3n -k4,4n > ${OUTPUTFILE4}
sed -i '1iGeneID\tChr\tStart\tEnd\tStrand' ${OUTPUTFILE4}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE4}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
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