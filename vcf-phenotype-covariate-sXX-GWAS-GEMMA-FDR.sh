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
##Filter samples in a VCF file (.vcf.gz) by minor allele frequency and missing data and create relatedness matrix using GEMMA.

##Input $1: Output location.
##Input $2: VCF file (.vcf.gz).
##Input $3: Tab-separated phenotype matrix (.txt; without header or row names), with samples as rows and phenotypes as columns.
##Input $4: Tab-separated covariate matrix (.txt; without header or row names), with samples as rows and covariates as columns.
##Input $5: Tab-separated relatedness matrix (.sXX.txt; without header or row names).
##Input $6: Comma-separated list of phenotype indices.
##Output: GWAS output file (.assoc.txt) and log file (.log.txt).

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
##    -J vcf-phenotype-covariate-sXX-GWAS-GEMMA-FDR-${<INPUTFILE>##*/} \
##    vcf-phenotype-covariate-sXX-GWAS-GEMMA-FDR.sh <OUTPUT LOCATION> <INPUT VCF FILE> <INPUT PHENOTYPE FILE> <INPUT COVARIATE FILE> <INPUT RELATEDNESS FILE> <LIST OF PHENOTYPE INDICES>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-phenotype-covariate-sXX-GWAS-GEMMA-FDR.sh") 

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
module load GEMMA
module load PDC/24.11
module load R/4.4.2-cpeGNU-24.11

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTVCFFILE=$(readlink -f $2)
INPUTPHENOTYPEFILE=$(readlink -f $3)
INPUTCOVARIATEFILE=$(readlink -f $4)
INPUTRELATEDNESSFILE=$(readlink -f $5)
PHENOTYPEINDICES=$6

if [[ -z "${PHENOTYPEINDICES}" ]] ; then
	NPHENOTYPES=$(head -n1 ${INPUTPHENOTYPEFILE} | tr '\t' '\n' | wc -l)
	PHENOTYPEINDICES=$(seq 1 ${NPHENOTYPES} | paste -sd,)
fi

############################################################################
##OUTPUT:

INPUTFILEPREFIX=$(echo "${INPUTVCFFILE##*/}" | sed 's/\.vcf$//' | sed 's/\.vcf\.gz$//' | sed 's/-job[0-9].*$//')

##Create output directory.
OUTPUTDIR=$(echo "${OUTPUTLOCATION}/${INPUTFILEPREFIX}.gemmagwas-job${JOBID}")
mkdir -p ${OUTPUTDIR}

OUTPUTLIST=$(echo "${OUTPUTDIR}/files.gemmagwas-job${JOBID}.txt")
ERRORFILE=$(echo "${OUTPUTDIR}/error.gemmagwas-job${JOBID}.txt")

TMPFILEPREFIX=$(echo "${OUTPUTDIR}/${INPUTFILEPREFIX}.gemma_plinkformat-job${JOBID}")

############################################################################
##ACTIONS:

##Filter VCF file and convert it to PLINK format.
#plink --vcf ${INPUTVCFFILE} --maf 0.05 --geno 0.05 --allow-extra-chr --chr-set 39 --make-bed --out ${TMPFILEPREFIX}
plink --vcf ${INPUTVCFFILE} --allow-extra-chr --chr-set 39 --make-bed --out ${TMPFILEPREFIX}

##Run GEMMA GWAS.
echo ${PHENOTYPEINDICES} | tr ',' '\n' | while read PHENOTYPEi ; do
	echo "Running GEMMA on phenotype ${PHENOTYPEi}."
	
	##Run association analysis.
	OUTPUTFILEPREFIX=$(echo "${INPUTFILEPREFIX}.gemma_gwasoutput_phenotype${PHENOTYPEi}-job${JOBID}")
	gemma \
		-bfile ${TMPFILEPREFIX} \
		-p ${INPUTPHENOTYPEFILE} \
		-n ${PHENOTYPEi} \
		-c ${INPUTCOVARIATEFILE} \
		-k ${INPUTRELATEDNESSFILE} \
		-lmm 4 \
		-outdir ${OUTPUTDIR} \
		-o ${OUTPUTFILEPREFIX}

	OUTPUTFILE="${OUTPUTDIR}/${OUTPUTFILEPREFIX}.assoc.txt"

	#Run R to perform FDR correction on PvalueLRT.
	TMPFILEFDR=$(echo "${OUTPUTDIR}/tmp.${PHENOTYPEi}.FDR-job${JOBID}.txt")
	paste <(cat ${OUTPUTFILE}) <(Rscript --vanilla -e "
	#Read input data.
	DATA<-read.table('${OUTPUTFILE}',header=TRUE,sep='\t',stringsAsFactors=FALSE)

	#Extract p-values.
	PVALUES<-as.numeric(DATA[,14])

	#Perform FDR correction.
	FDRPVALUES<-p.adjust(PVALUES,method='fdr')

	#Print FDR p-values.
	cat(FDRPVALUES,sep='\n')
	" | sed '1iPvalueFDR') > ${TMPFILEFDR}
	mv ${TMPFILEFDR} ${OUTPUTFILE}
	rm -f ${TMPFILEFDR}

	##Check for errors.
	if [[ ! -s "${OUTPUTFILE}" ]] ; then
		echo -e "ERROR:\t${OUTPUTDIR}/${OUTPUTFILEPREFIX}.assoc.txt"
		echo -e "${OUTPUTFILE}" >> ${ERRORFILE}
		echo -e "${OUTPUTFILE}\tERROR" >> ${OUTPUTLIST}
	else
		echo -e "${OUTPUTFILE}\tCLEAR" >> ${OUTPUTLIST}
	fi

done

##Remove temporary files.
rm -f ${TMPFILEPREFIX}*

############################################################################
##SAVE CONTROL FILES:

if [[ ! -s "${ERRORFILE}" ]]; then
	
	##README LOG ENTRY:
printf '%b\n' "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTVCFFILE}
Input file: ${INPUTPHENOTYPEFILE}
Input file: ${INPUTCOVARIATEFILE}
Input file: ${INPUTRELATEDNESSFILE}
Input phenotype indices: ${PHENOTYPEINDICES}
Output directory: ${OUTPUTDIR}
$(readlink -f $(find ${OUTPUTDIR}/*-job${JOBID}* -maxdepth 1) | sed 's/^/Output file: /g')
" >> $(echo "${OUTPUTLOCATION}/README.txt") 

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