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
##Filter samples in a VCF file (.vcf.gz) by minor allele frequency and missing data and create relatedness matrix using GEMMA.

##Input $1: Output location.
##Input $2: VCF file (.vcf.gz).
##Output: Standardized relatedness matrix file (.sXX.txt), log file (.log.txt) and Standardized relatedness matrix file with header (.sXXwithheader.txt).

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
##    -J vcf-relatednessmatrix-GEMMA-${<INPUTFILE>##*/} \
##    vcf-relatednessmatrix-GEMMA.sh <OUTPUT LOCATION> <INPUT VCF FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-relatednessmatrix-GEMMA.sh") 

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
module load plink
module load GEMMA

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)

############################################################################
##OUTPUT:

INPUTFILEPREFIX=$(echo "${INPUTFILE##*/}" | sed 's/\.vcf$//' | sed 's/\.vcf\.gz$//' | sed 's/-job[0-9].*$//')
OUTPUTFILEPREFIX=$(echo "${INPUTFILEPREFIX}.gemmarelatedness-job${JOBID}")
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.sXXwithheader.txt")

TMPPLINKPREFIX=$(echo "${OUTPUTLOCATION}/tmp1-job${JOBID}")
TMPPHENOTYPES=$(echo "${OUTPUTLOCATION}/tmp2-job${JOBID}.txt")
TMPHEADER=$(echo "${OUTPUTLOCATION}/tmp3-job${JOBID}.txt")
TMPROWS=$(echo "${OUTPUTLOCATION}/tmp4-job${JOBID}.txt")

############################################################################
##ACTIONS:

##Count chromosomes in VCF file.
CHRSET=$(bcftools index -s ${INPUTFILE} | cut -f 1 | wc -l)
##Filter VCF file and convert it to PLINK format.
plink --vcf ${INPUTFILE} --maf 0.05 --geno 0.05 --allow-extra-chr --chr-set ${CHRSET} --make-bed --out ${TMPPLINKPREFIX}

##Create dummy phenotype.
NSAMPLES=$(bcftools query -l ${INPUTFILE} | wc -l)
printf '1\n%.0s' $(seq ${NSAMPLES}) > ${TMPPHENOTYPES}

##Get relatedness matrix.
gemma \
	-bfile ${TMPPLINKPREFIX} \
	-p ${TMPPHENOTYPES} \
	-gk 2 \
	-outdir ${OUTPUTLOCATION} \
	-o ${OUTPUTFILEPREFIX}

##Create relatedness file with header (sample names).
bcftools query -l ${INPUTFILE} | tr '\n' '\t' > ${TMPHEADER}
echo "" >> ${TMPHEADER}
cat "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.sXX.txt" >> ${TMPHEADER}
bcftools query -l ${INPUTFILE} | sed '1iSAMPLE' > ${TMPROWS}
paste -d"\t" ${TMPROWS} ${TMPHEADER} > ${OUTPUTFILE}

##Remove temporary files.
rm -f $(readlink -f "${OUTPUTLOCATION}/tmp*-job${JOBID}*")

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.sXX.txt" ]] ; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Output file: ${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.sXX.txt
Output file: ${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.log.txt
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