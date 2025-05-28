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
##Phase VCF file (.vcf or .vcf.gz) using SHAPEIT5.

##Input $1: Output location.
##Input $2: VCF file (.vcf or .vcf.gz).
##Input $3: File containing a list of names of haploid chromosomes.
##Output: Phased VCF file (.vcf.gz) and its index file (.tbi).

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
##    -J vcf-phase-shapeit5-${<INPUTFILE>##*/} \
##    vcf-phase-shapeit5.sh <OUTPUT LOCATION> <INPUT VCF FILE> <HAPLOID CHR LIST FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-phase-shapeit5.sh") 

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
module load SHAPEIT/v5.1.1

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)
HAPLOIDLIST=$(readlink -f $3) ##List of haploid chromosomes.

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/\.bcf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILE1NAME=$(echo "${OUTPUTFILEPREFIX}.phase-job${JOBID}.vcf.gz") 
OUTPUTFILE2NAME=$(echo "${OUTPUTFILE1NAME}.tbi")
OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE1NAME}") 
OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE2NAME}") 

############################################################################
##ACTIONS:

##Create temporary output directory.
TMPDIR=$(echo "${OUTPUTLOCATION}/tmp-vcf-phase-shapeit5-job${JOBID}")
FILELIST=$(echo "${TMPDIR}/list-bcffiles.txt") 
mkdir -p ${TMPDIR}

##Phasing.
bcftools index -s ${INPUTFILE} | cut -f 1 | while read C ; do 
	echo "Processing ${C}."
	OUTPUTFILEXNAME=$(echo "${C}.bcf") 
	OUTPUTFILEX=$(echo "${TMPDIR}/${OUTPUTFILEXNAME}") 
	phase_common --input ${INPUTFILE} --region ${C} --filter-snp --haploids ${HAPLOIDLIST} --output-format bcf --output ${OUTPUTFILEX} --thread 10
	echo "${OUTPUTFILEX}" >> ${FILELIST}
done
sleep 5s

##Concatenate phased BCF files.
bcftools concat -O z --file-list ${FILELIST} --write-index=tbi -o ${OUTPUTFILE1}
sleep 5s

##Delete temporary output directory.
rm -rf ${TMPDIR}

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