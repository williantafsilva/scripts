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
##Create global consensus (IUPAC code) and individual consensus (alternative alleles in heterozygotes) sequences from a VCF file and the corresponding reference FASTA file.

##Input $1: Output location.
##Input $2: VCF file (.vcf or .vcf.gz).
##Input $3: Indexed reference genome FASTA file (.fa).
##Input $4: Target sequence range (chr:position-position).
##Input $5: Comma-separated string of sample names.
##Input $6: File tag (to be included in the output file name).
##Output: Global and individual sample consensus FASTA files (.fa).

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
##    -J vcf-consensusseq-${<INPUTFILE>##*/} \
##    vcf-consensusseq.sh <OUTPUT LOCATION> <INPUT VCF FILE> <REFERENCE FASTA FILE> <SEQUENCE RANGE> <SAMPLES> <FILE TAG>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-consensusseq.sh") 

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
REFFILE=$(readlink -f $3) ##Reference genome FASTA file.
SEQRANGE=$4 #Sequence range.
SAMPLELIST=$5 #List of sample names.
FILETAG=$6 #File tag.

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/\.bcf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILE1NAME=$(echo "${OUTPUTFILEPREFIX}.consensus${FILETAG}-global-job${JOBID}.fa") 
OUTPUTFILE2NAME=$(echo "${OUTPUTFILEPREFIX}.consensus${FILETAG}-individual-job${JOBID}.fa") 
OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE1NAME}") 
OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE2NAME}") 

############################################################################
##ACTIONS:

cd ${OUTPUTLOCATION}

#Global consensus.
samtools faidx ${REFFILE} ${SEQRANGE} | bcftools consensus --iupac-codes --samples ${SAMPLELIST} ${INPUTFILE} -p 'Consensus_global_IUPAC_code | ' --output ${OUTPUTFILE1}

#Individual consensus.
IFS=$','
for IND in ${SAMPLELIST} ; do
	samtools faidx ${REFFILE} ${SEQRANGE} | bcftools consensus --haplotype A --samples ${IND} ${INPUTFILE} -p "Consensus_altallele_${IND} | " >> ${OUTPUTFILE2} 
done
IFS=$'\n'

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE1}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Input reference file: ${REFFILE}
Input target sequence: ${SEQRANGE}
Input sample list: ${SAMPLELIST}
Input file tag: ${FILETAG}
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