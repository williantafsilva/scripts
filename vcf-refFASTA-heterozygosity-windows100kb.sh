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
##Calculate heterozygosity in 100kb non-overlapping windows, using the reference FASTA file to calculate window ranges.

##Input $1: Output location.
##Input $2: Reference FASTA file (.fa, .fasta, .fna).
##Input $3: VCF file (.vcf or .vcf.gz).
##Output: Text file (.txt) with columns: GenomicWindow, Chromosome, Start, End, Sample, nRefHom, nNonRefHom, nHets, HetProp.

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
##    -J vcf-refFASTA-heterozygosity-windows100kb-${<INPUTFILE>##*/} \
##    vcf-refFASTA-heterozygosity-windows100kb.sh <OUTPUT LOCATION> <INPUT REFERENCE FASTA FILE> <INPUT VCF FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-refFASTA-heterozygosity-windows100kb.sh") 

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
module load samtools
module load bedtools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFASTAFILE=$(readlink -f $2)
INPUTVCFFILE=$(readlink -f $3)

INPUTFILELOCATION=${INPUTVCFFILE%/*}
INPUTFILENAME=${INPUTVCFFILE##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/\.bcf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILENAME=$(echo "${OUTPUTFILEPREFIX}.het100kb-job${JOBID}.txt") 
OUTPUTFILE=$(echo "${OUTPUTLOCATION}/${OUTPUTFILENAME}") 

#Temporary files.
TMP1="${OUTPUTLOCATION}/tmp-het100kb-RefFastaIndex-job${JOBID}.fai"
TMP2="${OUTPUTLOCATION}/tmp-het100kb-RefFastaWindows-job${JOBID}.bed"
TMP3="${OUTPUTLOCATION}/tmp-het100kb-stats-job${JOBID}.txt"

############################################################################
##ACTIONS:

#Index reference genome FASTA file. 
samtools faidx --fai-idx ${TMP1} ${INPUTFASTAFILE}

#Create non-overlapping windows.
bedtools makewindows -g ${TMP1} -w 100000 > ${TMP2}

#Calculate statistics (including heterozygosity).
echo -e "GenomicWindow\tChromosome\tStart\tEnd\tSample\tnRefHom\tnNonRefHom\tnHets\tHetProp" > ${OUTPUTFILE}
cat ${TMP2} | while read L ; do
	WINDOWCHR="$(echo ${L} | cut -f1)"
	WINDOWSTART="$(echo ${L} | cut -f2)"
	WINDOWEND="$(echo ${L} | cut -f3)"
	GENOMICWINDOW="${WINDOWCHR}:${WINDOWSTART}-${WINDOWEND}"
	bcftools stats --samples - --regions ${GENOMICWINDOW} ${VCF} > ${TMP3}
	grep "^PSC" ${TMP3} | \
	cut -f3-6 | \
	awk -v OFS='\t' -v W="${GENOMICWINDOW}" -v WCHR="${WINDOWCHR}" -v WSTART="${WINDOWSTART}" -v WEND="${WINDOWEND}" '{ 
		if ($2 + $3 + $4 != 0) print W, WCHR, WSTART, WEND, $0, $4 / ( $2 + $3 + $4 ); else print W, WCHR, WSTART, WEND, $0, "NaN" 
	}' >> ${OUTPUTFILE}
done

rm -f ${TMP1} ${TMP2} ${TMP3}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input reference FASTA file: ${INPUTFASTAFILE}
Input VCF file: ${INPUTVCFFILE}
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