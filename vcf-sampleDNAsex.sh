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
##Change sample names in a VCF file (.vcf or .vcf.gz).

##Input $1: Output location.
##Input $2: VCF file (.vcf or .vcf.gz).
##Input $3: Sex-determining chromosome (Y in humans, W in chicken).
##Input $4: Essential sex chromosome (X in humans, Z in chicken).

##Output: File (.txt) with list of sample names, proportion of homozygotic sites in sex chromosome 1, 
##proportion of homozygotic sites in sex chromosome 2 and predicted genetic sex.

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
##    -J vcf-sampleDNAsex-${<INPUTFILE>##*/} \
##    vcf-sampleDNAsex.sh <OUTPUT LOCATION> <INPUT VCF FILE> <SEX CHR 1 NAME> <SEX CHR 2 NAME>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-sampleDNAsex.sh") 

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
SEXCHR1=$3
SEXCHR2=$4

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/\.bcf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILENAME=$(echo "${OUTPUTFILEPREFIX}.sampleDNAsex-job${JOBID}.txt") 
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}") 

############################################################################
##ACTIONS:

echo -e "Sample\t${SEXCHR1}(0/0)\t${SEXCHR1}(0/1)\t${SEXCHR2}(0/0)\t${SEXCHR2}(0/1)\tGeneticSex" > "${OUTPUTFILE}"

bcftools query -l ${INPUTFILE} | while read SAMPLENAME ; do

	NHOMOSEXCHR1=$(bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n]' -r ${SEXCHR1} -s ${SAMPLENAME} ${INPUTFILE} | grep -E "0/0|1/1|2/2|0\|0|1\|1|2\|2" | wc -l)
	NHETESEXCHR1=$(bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n]' -r ${SEXCHR1} -s ${SAMPLENAME} ${INPUTFILE} | grep -E "0/1|1/0|0/2|2/0|1/2|2/1|0\|1|1\|0|0\|2|2\|0|1\|2|2\|1" | wc -l)
	NHOMOSEXCHR2=$(bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n]' -r ${SEXCHR2} -s ${SAMPLENAME} ${INPUTFILE} | grep -E "0/0|1/1|2/2|0\|0|1\|1|2\|2" | wc -l)
	NHETESEXCHR2=$(bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n]' -r ${SEXCHR2} -s ${SAMPLENAME} ${INPUTFILE} | grep -E "0/1|1/0|0/2|2/0|1/2|2/1|0\|1|1\|0|0\|2|2\|0|1\|2|2\|1" | wc -l)
	if [[ "${NHOMOSEXCHR1}" -gt 0 || "${NHETESEXCHR1}" -gt 0 ]] ; then
		SAMPLESEX="Heterogametic"
	else
		SAMPLESEX="Homogametic"
	fi
	echo -e "${SAMPLENAME}\t${NHOMOSEXCHR1}\t${NHETESEXCHR1}\t${NHOMOSEXCHR2}\t${NHETESEXCHR2}\t${SAMPLESEX}" >> "${OUTPUTFILE}"

done

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Input sex-determining chromosome: ${SEXCHR1}
Input essential sex chromosome: ${SEXCHR2}
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