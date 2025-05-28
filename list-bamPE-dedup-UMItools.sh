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
##Deduplicate BAM files using UMI_tools.

##Input $1: Output location.
##Input $2: List (.txt) with paths to sorted and indexed BAM files aligned to coordinates (.bam).
##Output: Deduplicated BAM file (.bam), index file (.bai), log file (.log) and stats files (.tsv).

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
##    -J list-bamPE-dedup-UMItools-${<INPUTFILE>##*/} \
##    list-bamPE-dedup-UMItools.sh <OUTPUT LOCATION> <INPUT LIST FILE> 

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "list-bamPE-dedup-UMItools.sh") 

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
module load samtools
module load umi_tools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTLIST=$(readlink -f $2)

INPUTLISTLOCATION=${INPUTLIST%/*}
INPUTLISTNAME=${INPUTLIST##*/}

############################################################################
##OUTPUT:

OUTPUTLIST=$(echo "${OUTPUTLOCATION}/files.UMIdedup-job${JOBID}.txt")
ERRORFILE=$(echo "${OUTPUTLOCATION}/error.UMIdedup-job${JOBID}.txt")

############################################################################
##ACTIONS:

echo -e "INPUTFILE\tOUTPUTDIR\tOUTPUTFILE1\tOUTPUTFILE2\tERROR?" > ${OUTPUTLIST}
cat ${INPUTLIST} | while read INPUTFILE ; do
	
	#######################################
	##INPUT:

	INPUTFILEPREFIX=$(echo ${INPUTFILE##*/} | sed 's/\.bam$//' | sed 's/-job[0-9].*$//')

	#######################################
	##OUTPUT:

	OUTPUTFILEPREFIX=$(echo "${INPUTFILEPREFIX}.UMIdedup-job${JOBID}")

	##Create sample output directory.
	OUTPUTDIR=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}") 
	mkdir -p ${OUTPUTDIR}

	OUTPUTFILE1NAME=$(echo "${OUTPUTFILEPREFIX}.bam") 
	OUTPUTFILE2NAME=$(echo "${OUTPUTFILEPREFIX}.log") 
	OUTPUTFILE1=$(echo "${OUTPUTDIR}/${OUTPUTFILE1NAME}")
	OUTPUTFILE2=$(echo "${OUTPUTDIR}/${OUTPUTFILE2NAME}")
	OUTPUTFILEXPREFIX=$(echo "${OUTPUTDIR}/${OUTPUTFILEPREFIX}")

	#######################################
	##ACTIONS:
	
	echo "-----> Procesing file: ${INPUTFILE}"

	umi_tools dedup \
	-I ${INPUTFILE} \
	-S ${OUTPUTFILE1} \
	-L ${OUTPUTFILE2} \
	--paired \
	--buffer-whole-contig \
	--multimapping-detection-method=NH \
	--output-stats=${OUTPUTFILEXPREFIX} \
	--umi-separator=':'

	##Sort BAM file.
	samtools sort -o ${OUTPUTFILE1} ${OUTPUTFILE1}

	##Index BAM file.
	samtools index -b ${OUTPUTFILE1}

	##Check for errors.
	if [[ ! -s "${OUTPUTFILE1}" ]] || [[ ! -s "${OUTPUTFILE2}" ]] ; then
		echo -e "ERROR:\t${INPUTFILE}"
		echo -e "${INPUTFILE}" >> ${ERRORFILE}
		echo -e "${INPUTFILE}\t${OUTPUTDIR}\t${OUTPUTFILE1}\t${OUTPUTFILE2}\tERROR" >> ${OUTPUTLIST}
	else
		echo -e "${INPUTFILE}\t${OUTPUTDIR}\t${OUTPUTFILE1}\t${OUTPUTFILE2}\tCLEAR" >> ${OUTPUTLIST}
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
Output file: ${OUTPUTLIST}
$(cat ${OUTPUTLIST} | grep -v "^INPUTFILE" | cut -f2 | sed "s/^/Output directory: /g")
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

cat ${OUTPUTLIST} | grep -v "^INPUTFILE" | grep -vP "\tERROR" | while read L ; do
	
	OUTPUTDIR=$(echo ${L} | cut -f2)

	##README LOG ENTRY:
printf '%b\n' "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input list: ${INPUTLIST}
Output file: ${OUTPUTLIST}
Output directory: $(echo -e ${L} | cut -f2 | sed "s/^/Output directory: /g")
$(echo -e ${L} | cut -f3-4 | tr '\t' '\n' | sed "s/^/Output file: /g")
" >> $(echo "${OUTPUTDIR}/README.txt") 

	##COPY OF SCRIPT:
	cat $0 > $(echo "${OUTPUTDIR}/job${JOBID}.sh")
	##COPY OF SLURM FILE:
	cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTDIR}/job${JOBID}.out") 

done

exit 0