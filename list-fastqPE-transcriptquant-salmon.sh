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
##Align and index paired-end FastQ files to a reference genome using STAR.

##Input $1: Output location.
##Input $2: List (.txt) with three tab-separated columns with sample names and FASTQ files (R1\tR2). 
##Input $3: Directory to reference transcriptome Salmon index.
##Input $4: GTF file with genome annotations.
##Output: Directories with quantification files (quant.sf), log files and metadata files for each sample.

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
##    -J list-fastqPE-transcriptquant-salmon-${<INPUTFILE>##*/} \
##    list-fastqPE-transcriptquant-salmon.sh <OUTPUT LOCATION> <INPUT FILE> <REFERENCE TRANSCRIPTOME SALMON INDEX DIRECTORY>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "list-fastqPE-transcriptquant-salmon.sh") 

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
module load salmon

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTLIST=$(readlink -f $2)
SALMONINDEX=$(readlink -f $3)
GTFFILE=$(readlink -f $4)

INPUTLISTLOCATION=${INPUTLIST%/*}
INPUTLISTNAME=${INPUTLIST##*/}

############################################################################
##OUTPUT:

OUTPUTLIST=$(echo "${OUTPUTLOCATION}/files.salmonquant-job${JOBID}.txt")
ERRORFILE=$(echo "${OUTPUTLOCATION}/error.salmonquant-job${JOBID}.txt")
TRANSCRIPTTOGENEMAP=$(echo "${OUTPUTLOCATION}/transcripttogenemap.salmonquant-job${JOBID}.txt")

############################################################################
##ACTIONS:

##Create Salmon conversion file mapping transcript ID to gene ID from GTF file.
GTFTYPE=$(cat ${GTFFILE} | grep -v "#" | awk '$3=="transcript"' | head -n 1 | cut -f9 | tr -s ";" " " | awk '{print$3}' | sort | uniq | sed 's/"//g')
if [[ ${GTFTYPE} == "transcript_id" ]]; then
	echo "Input GTF file is a Gencode GTF."
	cat ${GTFFILE} | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$4"\t"$2}' | sort | uniq | sed 's/"//g' > ${TRANSCRIPTTOGENEMAP}
elif [[ ${GTFTYPE} == "gene_version" ]]; then
	echo "Input GTF file is a non-Gencode GTF (e.g., Ensembl)."
	cat ${GTFFILE} | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$6 "." $8"\t"$2 "." $4}' | sort | uniq | sed 's/"//g' > ${TRANSCRIPTTOGENEMAP}
fi

echo -e "INPUTFILE1\tINPUTFILE2\tOUTPUTDIR\tERROR?" > ${OUTPUTLIST}
cat ${INPUTLIST} | while read LINE ; do
	
	#######################################
	##INPUT:

	INPUTFILE1=$(echo ${LINE} | cut -f1)
	INPUTFILE2=$(echo ${LINE} | cut -f2)

	#######################################
	##OUTPUT:

	INPUTFILE1PREFIX=$(echo ${INPUTFILE1##*/} | sed 's/\.fastq.gz$//' | sed 's/\.fq.gz$//' | sed 's/-job[0-9].*$//')
	INPUTFILE2PREFIX=$(echo ${INPUTFILE2##*/} | sed 's/\.fastq.gz$//' | sed 's/\.fq.gz$//' | sed 's/-job[0-9].*$//')

	##Create consensus file name prefix.
	CONSENSUSFILEPREFIX=""
	while read C ; do
	    C1=$(echo ${C} | cut -f1)
	    C2=$(echo ${C} | cut -f2)
	    if [[ "${C1}" == "${C2}" ]] ; then
	        CONSENSUSFILEPREFIX+="${C1}"
	    else
	        CONSENSUSFILEPREFIX+="X"
	    fi
	done <<< "$(paste -d'\t' <(echo ${INPUTFILE1PREFIX} | grep -o .) <(echo ${INPUTFILE2PREFIX} | grep -o .))"

	##Create sample output directory.
	OUTPUTDIR=$(echo "${OUTPUTLOCATION}/${CONSENSUSFILEPREFIX}.salmonquant-job${JOBID}") 

	#######################################
	##ACTIONS:

	echo "-----> Processing files: ${INPUTFILE1}, ${INPUTFILE2}"

	##Run Salmon quantification.
	salmon quant \
    	--index ${SALMONINDEX} \
    	--libType A \
    	-1 ${INPUTFILE1} \
    	-2 ${INPUTFILE2} \
    	-o ${OUTPUTDIR} \
    	--validateMappings \
    	--geneMap ${TRANSCRIPTTOGENEMAP} \
    	--seqBias \
    	--gcBias \
    	--posBias \
    	--writeUnmappedNames
    	
	##Check for errors.
	if [[ ! -s "${OUTPUTDIR}/quant.sf" ]] || [[ ! -s "${OUTPUTDIR}/quant.genes.sf" ]]; then
		echo -e "ERROR:\t${INPUTFILE1}\t${INPUTFILE2}"
		echo -e "${INPUTFILE1}\t${INPUTFILE2}" >> ${ERRORFILE}
		echo -e "${INPUTFILE1}\t${INPUTFILE2}\t${OUTPUTDIR}\tERROR" >> ${OUTPUTLIST}
	else
		echo -e "${INPUTFILE1}\t${INPUTFILE2}\t${OUTPUTDIR}\tCLEAR" >> ${OUTPUTLIST}
	fi

done

############################################################################
##SAVE CONTROL FILES:

if [[ ! -s "${ERRORFILE}" ]]; then
	
	##README LOG ENTRY:
printf '%b\n' "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input list: ${INPUTLIST}
Input Salmon index directory: ${SALMONINDEX}
Input GTF file: ${GTFFILE}
Output file: ${OUTPUTLIST}
Output file: ${TRANSCRIPTTOGENEMAP}
$(cat ${OUTPUTLIST} | grep -v "^INPUTFILE1" | cut -f3 | sed "s/^/Output directory: /g")
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

cat ${OUTPUTLIST} | grep -v "^INPUTFILE1" | grep -vP "\tERROR" | while read L ; do
	
	OUTPUTDIR=$(echo ${L} | cut -f3)

	##README LOG ENTRY:
printf '%b\n' "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input list: ${INPUTLIST}
Input index directory: ${SALMONINDEX}
Input GTF file: ${GTFFILE}
Output file: ${OUTPUTLIST}
Output file: ${TRANSCRIPTTOGENEMAP}
Output directory: $(echo -e ${L} | cut -f3)
" >> $(echo "${OUTPUTDIR}/README.txt") 

	##COPY OF SCRIPT:
	cat $0 > $(echo "${OUTPUTDIR}/job${JOBID}.sh")
	##COPY OF SLURM FILE:
	cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTDIR}/job${JOBID}.out") 

done

exit 0