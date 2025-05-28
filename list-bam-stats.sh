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
##Calculate genome-wide statistics and coverage of a BAM file.

##Input $1: Output location.
##Input $2: List (.txt) with paths to sorted and indexed BAM files (.bam).
##Output: Genome-wide statistics and coverage data (.txt).

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
##    -J list-bam-stats-${<INPUTFILE>##*/} \
##    list-bam-stats.sh <OUTPUT LOCATION> <INPUT LIST FILE>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "list-bam-stats.sh") 

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
module load bamtools
module load samtools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTLIST=$(readlink -f $2)

INPUTLISTLOCATION=${INPUTLIST%/*}
INPUTLISTNAME=${INPUTLIST##*/}

############################################################################
##OUTPUT:

OUTPUTLIST=$(echo "${OUTPUTLOCATION}/files.bamstats-job${JOBID}.txt")
ERRORFILE=$(echo "${OUTPUTLOCATION}/error.bamstats-job${JOBID}.txt")

############################################################################
##ACTIONS:

echo -e "INPUTFILE\tOUTPUTDIR\tOUTPUTFILE1\tOUTPUTFILE2\tOUTPUTFILE3\tOUTPUTFILE4\tOUTPUTFILE5\tERROR?" > ${OUTPUTLIST}
cat ${INPUTLIST} | while read INPUTFILE ; do
	
	#######################################
	##INPUT:

	INPUTFILEPREFIX=$(echo ${INPUTFILE##*/} | sed 's/\.bam$//' | sed 's/-job[0-9].*$//')

	#######################################
	##OUTPUT:

	OUTPUTFILEPREFIX=$(echo "${INPUTFILEPREFIX}.bamstats-job${JOBID}")
	
	##Create sample output directory.
	OUTPUTDIR=$(echo "${OUTPUTLOCATION}/${OUTPUTFILEPREFIX}") 
	mkdir -p ${OUTPUTDIR}

	OUTPUTFILE1NAME=$(echo "${INPUTFILEPREFIX}.bamstatsbamtools-job${JOBID}.txt")
	OUTPUTFILE2NAME=$(echo "${INPUTFILEPREFIX}.bamstatssamtools-job${JOBID}.txt")
	OUTPUTFILE3NAME=$(echo "${INPUTFILEPREFIX}.bamstatsindex-job${JOBID}.txt")
	OUTPUTFILE4NAME=$(echo "${INPUTFILEPREFIX}.bamstatsgenomecoverage-job${JOBID}.txt")
	OUTPUTFILE5NAME=$(echo "${INPUTFILEPREFIX}.bamstatsgenomecoverageplots-job${JOBID}.txt")
	OUTPUTFILE1=$(echo "${OUTPUTDIR}/${OUTPUTFILE1NAME}") 
	OUTPUTFILE2=$(echo "${OUTPUTDIR}/${OUTPUTFILE2NAME}") 
	OUTPUTFILE3=$(echo "${OUTPUTDIR}/${OUTPUTFILE3NAME}") 
	OUTPUTFILE4=$(echo "${OUTPUTDIR}/${OUTPUTFILE4NAME}") 
	OUTPUTFILE5=$(echo "${OUTPUTDIR}/${OUTPUTFILE5NAME}") 

	#######################################
	##ACTIONS:
	
	echo "-----> Procesing file: ${INPUTFILE}"

	##General statistics.
	bamtools stats -in ${INPUTFILE} > ${OUTPUTFILE1}
	samtools stats ${INPUTFILE} > ${OUTPUTFILE2}
	samtools idxstats ${INPUTFILE} > ${OUTPUTFILE3}
	samtools coverage ${INPUTFILE} >  ${OUTPUTFILE4}

	##Genome-wide coverage histogram.
	samtools coverage --ascii --plot-depth --n-bins 25 ${INPUTFILE} > ${OUTPUTFILE5}
	echo "
	Average genome-wide coverage (samtools depth): " >> ${OUTPUTFILE5}
	samtools depth -a ${INPUTFILE} | awk '{sum+=$3; sumsq+=$3*$3} END {print "Average coverage = ",sum/NR; print "Std. Dev. = ",sqrt(sumsq/NR - (sum/NR)**2)}' >> ${OUTPUTFILE5}

	##Check for errors.
	if [[ ! -s "${OUTPUTFILE1}" ]] || [[ ! -s "${OUTPUTFILE2}" ]] || [[ ! -s "${OUTPUTFILE3}" ]] || [[ ! -s "${OUTPUTFILE4}" ]] || [[ ! -s "${OUTPUTFILE5}" ]] ; then
		echo -e "ERROR:\t${INPUTFILE}"
		echo -e "${INPUTFILE}" >> ${ERRORFILE}
		echo -e "${INPUTFILE}\t${OUTPUTDIR}\t${OUTPUTFILE1}\t${OUTPUTFILE2}\t${OUTPUTFILE3}\t${OUTPUTFILE4}\t${OUTPUTFILE5}\tERROR" >> ${OUTPUTLIST}
	else
		echo -e "${INPUTFILE}\t${OUTPUTDIR}\t${OUTPUTFILE1}\t${OUTPUTFILE2}\t${OUTPUTFILE3}\t${OUTPUTFILE4}\t${OUTPUTFILE5}\tCLEAR" >> ${OUTPUTLIST}
	fi

done

############################################################################
##SAVE CONTROL FILES:

if [[ ! -s "${ERRORFILE}" ]]; then ?????
	
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
$(echo -e ${L} | cut -f3-7 | tr '\t' '\n' | sed "s/^/Output file: /g")
" >> $(echo "${OUTPUTDIR}/README.txt") 

	##COPY OF SCRIPT:
	cat $0 > $(echo "${OUTPUTDIR}/job${JOBID}.sh")
	##COPY OF SLURM FILE:
	cat $(echo "${PATHTOMYSLURM}/slurm-${JOBID}.out") > $(echo "${OUTPUTDIR}/job${JOBID}.out") 

done

exit 0