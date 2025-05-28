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
##Split indexed VCF file (.vcf or .vcf.gz) by chromosome/contig.

##Input $1: Output location.
##Input $2: Indexed VCF file (.vcf or .vcf.gz).
##Output: Directory with VCF files (.vcf.gz) by chromosome/contig.

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
##    -J vcf-split-${<INPUTFILE>##*/} \
##    vcf-split.sh <OUTPUT LOCATION> <INPUT VCF FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-split.sh") 

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

OUTPUTDIRPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/\.bcf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTDIRNAME=$(echo "${OUTPUTDIRPREFIX}.splitbychr-job${JOBID}") 
OUTPUTDIR=$(echo "${OUTPUTLOCATION}/${OUTPUTDIRNAME}") 

############################################################################
##ACTIONS:

##Create output directory.
mkdir -p ${OUTPUTDIR}

##Split VCF file.
while read C ; do 
    echo "Processing ${C}."
    OUTPUTFILEXNAME=$(echo "${C}.splitbychr-job${JOBID}.vcf.gz") 
    OUTPUTFILEX=$(echo "${OUTPUTDIR}/${OUTPUTFILEXNAME}") 
    bcftools view -O z -o ${OUTPUTFILEX} ${INPUTFILE} "${C}" 
    OUTPUTFILEXLIST=$(echo "${OUTPUTFILEXLIST}
Output file: ${OUTPUTFILEX}")
done <<< "$(bcftools index -s ${INPUTFILE} | cut -f 1)"

############################################################################
##SAVE CONTROL FILES:

if [[ $(ls ${OUTPUTDIR} | wc -l) -gt 0 ]]; then
    
    ##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Output directory: ${OUTPUTDIR}
" >> $(echo "${OUTPUTLOCATION}/README.txt") 

echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}${OUTPUTFILEXLIST}
" >> $(echo "${OUTPUTDIR}/README.txt") 

    ##COPY OF SCRIPT:
    cat $0 > $(echo "${OUTPUTLOCATION}/job${JOBID}.sh")
    cat $0 > $(echo "${OUTPUTDIR}/job${JOBID}.sh")
    ##COPY OF SLURM FILE:
    cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTLOCATION}/job${JOBID}.out") 
    cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTDIR}/job${JOBID}.out") 

else

    ##COPY OF SCRIPT:
    cat $0 > $(echo "${OUTPUTLOCATION}/job${JOBID}.sh")
    cat $0 > $(echo "${OUTPUTDIR}/job${JOBID}.sh")
    ##COPY OF SLURM FILE:
    cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTLOCATION}/job${JOBID}.err") 
    cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTDIR}/job${JOBID}.err")

fi

exit 0