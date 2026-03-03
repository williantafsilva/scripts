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
OUTPUTFILE1=$(echo "${OUTPUTFILEPREFIX}.bed")
OUTPUTFILE2=$(echo "${OUTPUTFILEPREFIX}.bin")
OUTPUTFILE3=$(echo "${OUTPUTFILEPREFIX}.fam")
OUTPUTFILE4=$(echo "${OUTPUTFILEPREFIX}.ld.gz")
OUTPUTFILE5=$(echo "${OUTPUTFILEPREFIX}.ldperdist.txt")
OUTPUTFILE6=$(echo "${OUTPUTFILEPREFIX}.ldper10kbdistbin.txt")
OUTPUTFILE7=$(echo "${OUTPUTFILEPREFIX}.meanldper10kbdistbin.txt")

############################################################################
##ACTIONS:

#Define chromosomes.
CHRNUMBER=$(bcftools index -s ${INPUTFILE} | cut -f1 | grep -E "^[0-9].*|^W|^Z|^X|^Chr[0-9].*|^chr[0-9].*|^ChrX|^ChrY|^ChrZ|^ChrW|^chrX|^chrY|^chrZ|^chrW" | wc -l)
CHRSET=$(bcftools index -s ${INPUTFILE} | cut -f1 | grep -E "^[0-9].*|^W|^Z|^X|^Chr[0-9].*|^chr[0-9].*|^ChrX|^ChrY|^ChrZ|^ChrW|^chrX|^chrY|^chrZ|^chrW" | paste -sd,)

echo "Compute pairwise LD (r^2)."

plink --vcf ${INPUTFILE} \
--double-id \
--chr ${CHRSET} \
--chr-set ${CHRNUMBER} \
--allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.05 \
--geno 0.1 \
--mind 0.5 \
--thin 0.5 \
--r2 gz \
--ld-window 100 \
--ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed \
--out ${OUTPUTFILEPREFIX}

echo "Calculate LD per chromosome and distance."

zcat ${OUTPUTFILE4} | \
awk '
BEGIN {
    print "Chromosome\tDistance_bp\tr2"
}
NR > 1 && $1 == $4 {
    dist = ($5 > $2 ? $5 - $2 : $2 - $5)
    print $1 "\t" dist "\t" $7
}' | sort -k1,1 -k2,2n > ${OUTPUTFILE5}

echo "Calculate LD per 10kb distance bin."

zcat ${OUTPUTFILE4} | \
awk '
BEGIN {
    print "DistanceBin_bp\tr2"
}
NR > 1 && $1 == $4 {
    dist = ($5 > $2 ? $5 - $2 : $2 - $5)
    bin = int(dist / 10000) * 10000
    print bin "\t" $7
}' | sort -n > ${OUTPUTFILE6}

echo "Calculate average LD per 10kb distance bin."

zcat ${OUTPUTFILE4} | \
awk 'NR > 1 && $1 == $4 {
    dist = ($5 > $2 ? $5 - $2 : $2 - $5)
    bin = int(dist / 10000) * 10000
    sum[bin] += $7
    count[bin]++
}
END {
    print "Distance_bp\tMean_r2\tN_pairs"
    for (b in sum)
        print b "\t" sum[b]/count[b] "\t" count[b]
}' | sort -n > ${OUTPUTFILE7}

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
Output file: ${OUTPUTFILE5}
Output file: ${OUTPUTFILE6}
Output file: ${OUTPUTFILE7}
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