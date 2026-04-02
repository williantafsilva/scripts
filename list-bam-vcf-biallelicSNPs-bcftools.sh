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
##Variant calling (only biallelic SNPs) with filtering using bcftools.

##Input $1: Output location.
##Input $2: Prefix of output VCF file.
##Input $2: Indexed reference genome FASTA file (*.fasta).
##Input $3: List (.txt) with sorted and indexed BAM files (*.bam).
##Output: Compressed VCF file (.vcf.gz).

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
##    -J list-bam-vcf-biallelicSNPs-bcftools-${<INPUTFILE>##*/} \
##    list-bam-vcf-biallelicSNPs-bcftools.sh <OUTPUT LOCATION> <PREFIX OF OUTPUT VCF FILE> <INPUT REFERENCE FASTA FILE> <INPUT BAM LIST>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "list-bam-vcf-biallelicSNPs-bcftools.sh") 

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
module load bcftools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
OUTPUTFILEPREFIX=$(readlink -f $2)
INPUTREFFASTA=$(readlink -f $3)
INPUTBAMLIST=$(readlink -f $4)

############################################################################
##OUTPUT:

OUTPUTFILE="${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.bcftoolsvcf-job${JOBID}.vcf.gz"

############################################################################
##ACTIONS:

#Index reference FASTA file.
samtools faidx ${INPUTREFFASTA}

#Variant calling.
bcftools mpileup \
	--fasta-ref=${INPUTREFFASTA} \
	--no-BAQ \
	--annotate DP,AD \
	--min-MQ 20 \
	--min-BQ 20 \
	--max-depth 5000 \
	--output-type u \
	--bam-list ${INPUTBAMLIST} | \
bcftools call \
	--multiallelic-caller \
	--variants-only \
	--output-type u | \
bcftools view \
	--types snps \
	--min-alleles 2 \
	--max-alleles 2 | \
bcftools filter \
	--soft-filter LowQual \
	--include 'QUAL>=30 && DP>=10 && MQ>=30' | \
bcftools view \
	--apply-filters PASS \
	--output-type z \
	--write-index=tbi \
	--output ${OUTPUTFILE}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input output file prefix: ${OUTPUTFILEPREFIX}
Input reference FASTA file: ${INPUTREFFASTA}
Input BAM list file: ${INPUTBAMLIST}
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