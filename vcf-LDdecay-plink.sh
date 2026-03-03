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
##Calculate average genome-wide linkage disequilibrium decay (LD decay) curve from a VCF file.

##Input $1: Output location.
##Input $2: Indexed VCF file (.vcf.gz).
##Output: Table containing columns: Distance_bp, Average_r2 

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
##    -J vcf-LDdecay-plink-${<INPUTFILE>##*/} \
##    vcf-LDdecay-plink.sh <OUTPUT LOCATION> <INPUT VCF FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-LDdecay-plink.sh")

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
BIN_SIZE=10000 #Bin size.

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

INPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILEPREFIX=$(echo "${OUTPUTLOCATION}/${INPUTFILEPREFIX}.plinkLDdecay-job${JOBID}") 
OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.bed")
OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.bin")
OUTPUTFILE3=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.fam")
OUTPUTFILE4=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}.ld.gz")

############################################################################
##ACTIONS:

#Count chromosomes.
CHRSET=$(bcftools index -s ${INPUTFILE} | cut -f 1 | wc -l)

echo "Compute pairwise LD (r^2)."

plink --vcf ${INPUTFILE} \
--double-id \
--chr-set ${CHRSET} \
--allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.05 --geno 0.1 --mind 0.5 \
--thin 0.1 \
--r2 gz \
--ld-window 100 \
--ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed \
--out ${OUTPUTFILEPREFIX}

echo "Bin LD by distance"

#awk -v BIN=${BIN_SIZE} '
#BEGIN {
#    print "distance_bp\tmean_r2"
#}
#NR > 1 {
#    dist = ($5 > $2 ? $5 - $2 : $2 - $5)
#    bin = int(dist / BIN) * BIN
#    sum[bin] += $7
#    count[bin]++
#}
#END {
#    for (b in sum) {
#        print b "\t" sum[b]/count[b]
#    }
#}
#' ${OUTPUTFILEPREFIX}_ld.ld | sort -n > ${OUTPUTFILEPREFIX}_ld_decay.txt

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