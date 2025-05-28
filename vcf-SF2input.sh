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
##Generate SweepFinder2 input files (.SF2input).

##Input $1: Output location.
##Input $2: VCF file (.vcf.gz) split per chromosome/contig.
##Output: SweepFinder2 input file (.SF2input).

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
##    -J vcf-SF2input-${<INPUTFILE>##*/} \
##    vcf-SF2input.sh <OUTPUT LOCATION> <INPUT VCF FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-SF2input.sh") 

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
module load vcftools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/\.bcf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILENAME=$(echo "${OUTPUTFILEPREFIX}.SF2format-job${JOBID}.SF2input") 
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}") 

############################################################################
##ACTIONS:

##Get allele counts from the VCF files.
vcftools --gzvcf ${INPUTFILE} --counts2 --out $(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}-job${JOBID}")
sleep 5s

##Format output files to SweepFinder2 standard (and remove monomorphic sites).
tail -n+2 $(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}-job${JOBID}.frq.count") | awk -v OFS="\t" '{print $2,$6,$4,"1"}' | sed '1 i\position\tx\tn\tfolded' > ${OUTPUTFILE} 
##tail -n+2 $(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}-job${JOBID}.frq.count") | awk -v OFS="\t" '{print $2,$6,$4,"1"}' | sed '1 i\position\tx\tn\tfolded' | awk -F'\t' '$2!=$3' > ${OUTPUTFILE} 
sleep 5s

##Remove .frq.count files.
rm -f $(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}-job${JOBID}.frq.count") 

##Remove .log files.
rm -f $(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}-job${JOBID}.log") 

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
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