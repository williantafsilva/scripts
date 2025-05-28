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
##Convert VCF file into genotype matrix file with genotypes encoded as 1, 2, 3, 4, ...

##Input $1: Output location.
##Input $2: Indexed VCF file (.vcf.gz).
##Output: Genotype matrix file with genotypes in VCF format (.txt), genotype matrix file with translated genotypes (.txt), and
##genotype matrix file with genotypes in numeric format (.txt).

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
##    -J vcf-genotypematrix-MatrixEQTL-${<INPUTFILE>##*/} \
##    vcf-genotypematrix-MatrixEQTL.sh <OUTPUT LOCATION> <INPUT VCF FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-genotypematrix-MatrixEQTL.sh")

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

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf\.gz$//' | sed 's/-job[0-9].*$//')
OUTPUTFILE1NAME=$(echo "${OUTPUTFILEPREFIX}.genotypematrixGTformat-job${JOBID}.txt") 
OUTPUTFILE2NAME=$(echo "${OUTPUTFILEPREFIX}.genotypematrixNUMformat-job${JOBID}.txt") 
OUTPUTFILE3NAME=$(echo "${OUTPUTFILEPREFIX}.genotypematrixTGTformat-job${JOBID}.txt") 
OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE1NAME}")
OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE2NAME}")
OUTPUTFILE3=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE3NAME}")

############################################################################
##ACTIONS:

##Create output files with headers.
{
    echo -ne "SNP\t$(bcftools query -l "${INPUTFILE}" | tr '\n' '\t' | sed 's/\t$/\n/')"
    echo ""
} > "${OUTPUTFILE1}"

{
    echo -ne "SNP\t$(bcftools query -l "${INPUTFILE}" | tr '\n' '\t' | sed 's/\t$/\n/')"
    echo ""
} > "${OUTPUTFILE3}"
bcftools query -f '%CHROM:%POS[\t%TGT]\n' "${INPUTFILE}" >> "${OUTPUTFILE3}"

{
    echo -ne "SNP\t$(bcftools query -l "${INPUTFILE}" | tr '\n' '\t' | sed 's/\t$/\n/')"
    echo ""
} > "${OUTPUTFILE2}"

##Extract genotype data and write to output files.
bcftools query -f '%CHROM:%POS[\t%GT]\n' "${INPUTFILE}" | \
awk -v OUTPUTGT="${OUTPUTFILE1}" -v OUTPUTNUM="${OUTPUTFILE2}" '
{
    #SNP ID.
    snp_id = $1;

    #Process each genotype column.
    gt_line = snp_id;
    num_line = snp_id;
    
    for (i=2; i<=NF; i++) {
        if ($i ~ /\.\/\./) { 
            gt_line = gt_line "\tNA"; 
            num_line = num_line "\tNA";
        } else {
            gt_line = gt_line "\t" $i;
            split($i, gt, /[\/|]/);
            num_line = num_line "\t" (gt[1] + gt[2]);
        }
    }

    #Append data to output files.
    print gt_line >> OUTPUTGT;
    print num_line >> OUTPUTNUM;
}'

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE2}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Output file: ${OUTPUTFILE1}
Output file: ${OUTPUTFILE2}
Output file: ${OUTPUTFILE3}
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