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
##Generate OmegaPlus output files.

##Input $1: Output location.
##Input $2: VCF file (.vcf.gz) split per chromosome/contig.
##Output: OmegaPlus_Info and OmegaPlus_Report files.

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
##    -J vcf-OmegaPlus-${<INPUTFILE>##*/} \
##    vcf-OmegaPlus.sh <OUTPUT LOCATION> <INPUT VCF FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-OmegaPlus.sh") 

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
module load htslib

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

INPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILESUFFIX=$(echo "${INPUTFILEPREFIX}.OmegaPlus-job${JOBID}.OmegaPlusoutput") 

############################################################################
##ACTIONS:

##Go to output location.
cd ${OUTPUTLOCATION}

##Decompress input file.
INPUTFILEX=$(echo ${INPUTFILENAME} | sed -E 's/(.*).vcf.gz/\1.vcf/')
bgzip -dc ${INPUTFILE} > ${INPUTFILEX}

##Calculate parameters.
SEQINTERVAL=10000 ##Genetic interval (bp).
POSMAX=$(bcftools query -f '%POS\n' ${INPUTFILE} | awk 'NR==1{max = $1 + 0; next} {if ($1 > max) max = $1;} END {print max}') ##Max position.
POSMIN=$(bcftools query -f '%POS\n' ${INPUTFILE} | awk 'NR==1{min = $1 + 0; next} {if ($1 < min) min = $1;} END {print min}') ##Min position.
NSITES=$(echo "$(((${POSMAX}-${POSMIN})/${SEQINTERVAL}+1))") ##Number of sites (grid size).
ALIGNLENGTH=$(echo "$((${POSMAX}-${POSMIN}))") ##Alignment length (default 200000).

##Run OmegaPlus.
OmegaPlus \
	-name ${OUTPUTFILESUFFIX} \
	-input ${INPUTFILEX} \
	-ld D \
	-grid ${NSITES} \
	-length ${ALIGNLENGTH} \
	-minwin 5000 \
	-maxwin 100000 \
	-seed 15243
sleep 1s

##Remove "//1" from the reportfile.
sed -i '/\/\//d' $(echo "OmegaPlus_Report.${OUTPUTFILESUFFIX}") 

##Remove empty lines from the report file.
sed -i '/^$/d' $(echo "OmegaPlus_Report.${OUTPUTFILESUFFIX}")

##Delete .vcf file.
rm -f ${INPUTFILEX}

############################################################################
##SAVE CONTROL FILES:

if [[ -s $(echo "OmegaPlus_Report.${OUTPUTFILESUFFIX}") ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Output file: OmegaPlus_Info.${OUTPUTFILESUFFIX}
Output file: OmegaPlus_Report.${OUTPUTFILESUFFIX}
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