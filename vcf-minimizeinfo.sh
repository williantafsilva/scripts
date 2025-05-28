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
##Minimize VCF file (.vcf or .vcf.gz), removing non-essential header information and keeping only the genotypes of each sample.

##Input $1: Output location.
##Input $2: VCF file (.vcf or .vcf.gz).
##Output: VCF file (.vcf.gz) and index file (.tbi).

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
##    -J vcf-minimizeinfo-${<INPUTFILE>##*/} \
##    vcf-minimizeinfo.sh <OUTPUT LOCATION> <INPUT VCF FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-minimizeinfo.sh") 

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

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/\.bcf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILE1NAME=$(echo "${OUTPUTFILEPREFIX}.minimize-job${JOBID}.vcf.gz") 
OUTPUTFILE2NAME=$(echo "${OUTPUTFILE1NAME}.tbi")
OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE1NAME}") 
OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE2NAME}") 

TMPFILE1=$(echo "${OUTPUTLOCATION}/tmp.header-job${JOBID}.txt")
TMPFILE2=$(echo "${OUTPUTLOCATION}/tmp.data-job${JOBID}.vcf.gz")

############################################################################
##ACTIONS:

##Get essential header information.
bcftools view -h "${INPUTFILE}" | grep "^##fileformat\|^##FILTER=<ID=PASS\|^##INFO=<ID=AC\|^##INFO=<ID=AN\|^##FORMAT=<ID=GT\|^##contig=<ID=" > "${TMPFILE1}"
echo "##fileDate=${RUNDATE}" >> "${TMPFILE1}"
echo "##source=bcftools annotate -x ^INFO/AC,^INFO/AN,^FORMAT/GT --output-type z --output ${OUTPUTFILE1} ${INPUTFILE}" >> "${TMPFILE1}"
bcftools view -h "${INPUTFILE}" | grep "^#CHROM" >> "${TMPFILE1}"

##Remove all INFO fields except for AC (allele count) and AN (allele number) and all FORMAT fields except for GT (genotype).
bcftools annotate -x ^INFO/AC,^INFO/AN,^FORMAT/GT --output-type z --output ${TMPFILE2} ${INPUTFILE} 
bcftools reheader --header ${TMPFILE1} --output ${OUTPUTFILE1} ${TMPFILE2}
bcftools index --tbi ${OUTPUTFILE1}

##Remove temporary files.
rm -f ${TMPFILE1} ${TMPFILE2}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE1}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Output file: ${OUTPUTFILE1}
Output file: ${OUTPUTFILE2}
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