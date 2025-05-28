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
##Align FASTA sequences using ClustalOmega.

##Input $1: Output location.
##Input $2: Indexed phased VCF file (.vcf.gz).
##Output: Files with PCA eigenvalues (.eigenval) and eigenvectors (.eigenvec).

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
##    -J vcf-pca-plink-${<INPUTFILE>##*/} \
##    vcf-pca-plink.sh <OUTPUT LOCATION> <INPUT VCF FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-pca-plink.sh")

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

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

INPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILEPREFIX=$(echo "${OUTPUTLOCATION}/${INPUTFILEPREFIX}.plinkpca-job${JOBID}") 
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${INPUTFILEPREFIX}.plinkpcacomplete-job${JOBID}.txt")

TMPFILEPREFIX=$(echo "tmp-job${JOBID}")
TMPFILEPREFIX=$(echo "${OUTPUTLOCATION}/${TMPFILEPREFIX}") 

############################################################################
##ACTIONS:

#Count chromosomes.
CHRSET=$(bcftools index -s ${INPUTFILE} | cut -f 1 | wc -l)

#Convert VCF to PLINK format.
plink --vcf ${INPUTFILE} --chr-set ${CHRSET} --allow-extra-chr --make-bed --out ${TMPFILEPREFIX}

#Perform LD pruning.
plink --bfile ${TMPFILEPREFIX} --chr-set ${CHRSET} --allow-extra-chr --indep-pairwise 50 5 0.2 --out "${TMPFILEPREFIX}_pruned"
plink --bfile ${TMPFILEPREFIX} --chr-set ${CHRSET} --allow-extra-chr --extract "${TMPFILEPREFIX}_pruned.prune.in" --make-bed --out "${TMPFILEPREFIX}_filtered"

#Perform PCA.
plink --bfile "${TMPFILEPREFIX}_filtered" --chr-set ${CHRSET} --allow-extra-chr --pca 20 --out "${OUTPUTFILEPREFIX}"

#Delete temporary files.
rm -f $(find ${TMPFILEPREFIX}*)

#Format main output file.
echo "SAMPLE $(printf '%s ' PC{1..20})" | tr ' ' '\t' > ${OUTPUTFILE}
cat "${OUTPUTFILEPREFIX}.eigenvec" | cut -d' ' -f2- | tr ' ' '\t' >> ${OUTPUTFILE}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILEPREFIX}.eigenvec" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Output file: ${OUTPUTFILEPREFIX}.eigenval
Output file: ${OUTPUTFILEPREFIX}.eigenvec
Output file: ${OUTPUTFILEPREFIX}.nosex
Output file: ${OUTPUTFILEPREFIX}.log
output file: ${OUTPUTFILE}
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