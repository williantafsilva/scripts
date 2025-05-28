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
##Output: BED files (.bed) for genes, transcripts, exons (with exon ID) and exons (with gene ID).

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
##    -J gtf-bed-${<INPUTFILE>##*/} \
##    gtf-bed.sh <OUTPUT LOCATION> <INPUT GTF FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "gtf-bed.sh") 

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
OUTPUTFILE1NAME=$(echo "${OUTPUTFILEPREFIX}.gtftobedgenes-job${JOBID}.bed") 
OUTPUTFILE2NAME=$(echo "${OUTPUTFILEPREFIX}.gtftobedtranscripts-job${JOBID}.bed") 
OUTPUTFILE3NAME=$(echo "${OUTPUTFILEPREFIX}.gtftobedexons-job${JOBID}.bed") 
OUTPUTFILE4NAME=$(echo "${OUTPUTFILEPREFIX}.gtftobedexonsgeneid-job${JOBID}.bed") 
OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE1NAME}") 
OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE2NAME}") 
OUTPUTFILE3=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE3NAME}") 
OUTPUTFILE4=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE4NAME}") 

############################################################################
##ACTIONS:

##Convert GTF file to BED file for genes.
paste -d' ' <(grep -P "\tgene\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$6,".",$4 }') \
<(grep -P "\tgene\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | sed 's/.*\(ENSGALG[0-9]*\t\).*/\1/g') | \
awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6}' | \
sort -V -k1,1 -k2,2n -k3,3n > ${OUTPUTFILE1}

##Convert GTF file to BED file for transcripts.
paste -d' ' <(grep -P "\ttranscript\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$10,".",$4 }') \
<(grep -P "\ttranscript\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | sed 's/.*\(ENSGALT[0-9]*\t\).*/\1/g') | \
awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6}' | \
sort -V -k1,1 -k2,2n -k3,3n > ${OUTPUTFILE2}

##Convert GTF file to BED file for exons (with exon ID).
paste -d' ' <(grep -P "\texon\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$6,".",$4 }') \
<(grep -P "\texon\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | sed 's/.*\(ENSGALE[0-9]*\t\).*/\1/g') | \
awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6}' | \
sort -V -k1,1 -k2,2n -k3,3n > ${OUTPUTFILE3}

##Convert GTF file to BED file for exons (with gene ID).
paste -d' ' <(grep -P "\texon\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$6,".",$4 }') \
<(grep -P "\texon\t" ${INPUTFILE} | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | sed 's/.*\(ENSGALG[0-9]*\t\).*/\1/g') | \
awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6}' | \
sort -V -k1,1 -k2,2n -k3,3n > ${OUTPUTFILE4}

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