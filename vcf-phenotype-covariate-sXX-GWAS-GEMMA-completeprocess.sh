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
##Input $3: Tab-separated phenotype/covariate matrix (.txt; WITH header and row names), with samples as rows and traits/covariates as columns.
##Input $4: Tab-separated relatedness matrix (.sXX.txt; WITH header and row names).
##Input $5: Comma-separated list of samples.
##Input $6: Comma-separated list of phenotype names.
##Input $7: Comma-separated list of covariate names.
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
##    -J vcf-phenotype-covariate-sXX-GWAS-GEMMA-${<INPUTFILE>##*/} \
##    vcf-phenotype-covariate-sXX-GWAS-GEMMA-completeprocess.sh <OUTPUT LOCATION> <INPUT VCF FILE> <INPUT PHENOTYPE FILE> <INPUT COVARIATE FILE> <INPUT RELATEDNESS FILE> <LIST OF PHENOTYPE INDICES>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-phenotype-covariate-sXX-GWAS-GEMMA-completeprocess.sh") 

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

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTVCFFILE=$(readlink -f $2)
INPUTPHENOTYPEFILE=$(readlink -f $3)
INPUTRELATEDNESSFILE=$(readlink -f $4)
LIST_SAMPLES=$5
LIST_TRAITS=$6
LIST_COVARIATES=$7

############################################################################
##OUTPUT:

INPUTFILEPREFIX=$(echo "${INPUTVCFFILE##*/}" | sed 's/\.vcf$//' | sed 's/\.vcf\.gz$//' | sed 's/-job[0-9].*$//')

##Create output directory.
OUTPUTDIR=$(echo "${OUTPUTLOCATION}/${INPUTFILEPREFIX}.gemmagwas-job${JOBID}")
mkdir -p ${OUTPUTDIR}

##Output files.
OUTPUTFILE1="${OUTPUTDIR}/genotypes.gemmagwas-job${JOBID}.vcf.gz"
OUTPUTFILE2="${OUTPUTDIR}/list-samples.gemmagwas-job${JOBID}.txt"
OUTPUTFILE3="${OUTPUTDIR}/list-targettraits.gemmagwas-job${JOBID}.txt"
OUTPUTFILE4="${OUTPUTDIR}/list-covariates.gemmagwas-job${JOBID}.txt"
OUTPUTFILE5="${OUTPUTDIR}/matrix-phenotypes.gemmagwas-job${JOBID}.txt"
OUTPUTFILE6="${OUTPUTDIR}/matrix-covariates.gemmagwas-job${JOBID}.txt"
OUTPUTFILE7="${OUTPUTDIR}/matrix-relatedness.gemmagwas-job${JOBID}.sXX.txt"

OUTPUTLIST=$(echo "${OUTPUTDIR}/files.gemmagwas-job${JOBID}.txt")
ERRORFILE=$(echo "${OUTPUTDIR}/error.gemmagwas-job${JOBID}.txt")

##Temporary files.
TMP1="${OUTPUTDIR}/tmp1-phenotypes-job${JOBID}.txt"
TMP2="${OUTPUTDIR}/tmp2-phenotypes-transpose-job${JOBID}.txt"
TMP3="${OUTPUTDIR}/tmp3-targettraits-job${JOBID}.txt"
TMP4="${OUTPUTDIR}/tmp4-covariates-job${JOBID}.txt"
TMP5="${OUTPUTDIR}/tmp5-relatedness-job${JOBID}.txt"
TMP6="${OUTPUTDIR}/tmp6-relatedness-job${JOBID}.txt"
TMP7="${OUTPUTDIR}/tmp7-relatedness-job${JOBID}.txt"
TMP_PLINKFILEPREFIX="${OUTPUTDIR}/${INPUTFILEPREFIX}.gemma_plinkformat-job${JOBID}"
TMP8="${OUTPUTDIR}/tmp8-PvalueFDR-job${JOBID}.txt"
TMP9="${OUTPUTDIR}/tmp9-PvalueFDR-job${JOBID}.txt"
TMPDIR1="${OUTPUTDIR}/tmpdir1-PvalueFDR-job${JOBID}"
TMPDIR2="${OUTPUTDIR}/tmpdir2-PvalueFDR-job${JOBID}"
mkdir -p ${TMPDIR1}
mkdir -p ${TMPDIR2}

############################################################################
##ACTIONS:

##Select samples in VCF file.
##Select variants based on minor allele frequency (MAF>0.5) and fraction of samples with missing genotype (F_MISSING<0.05).
echo "##################################################"
echo "Select samples."
echo "Select variants based on minor allele frequency (MAF>0.5) and fraction of samples with missing genotype (F_MISSING<0.05)."
echo "Remove non-chromosomal contigs."
CHRSET_GAL7B="$(seq 1 39 | tr '\n' ',')$(printf '%s,' Chr{1..39})$(printf '%s,' chr{1..39})W,Z,ChrW,ChrZ,chrW,chrZ"
bcftools view \
  --samples ${LIST_SAMPLES} \
  --regions ${CHRSET_GAL7B} \
  --min-ac=1 \
  --output-type u \
  ${INPUTVCFFILE} | \
  bcftools view \
  --include 'MAF>0.05 & F_MISSING<0.05' \
  --output-type z \
  --write-index=tbi \
  --output ${OUTPUTFILE1}

##Get list of samples with genotypes.
bcftools query -l ${OUTPUTFILE1} > ${OUTPUTFILE2}

##Define order of samples in the phenotype file.
echo "##################################################"
echo "Define order of samples in the phenotype file."
head -n1 ${INPUTPHENOTYPEFILE} > ${TMP1}
cat ${OUTPUTFILE2} | while read SAMPLE ; do
    grep -P "^${SAMPLE}\t" ${INPUTPHENOTYPEFILE} >> ${TMP1}
done

##Transpose phenotype file.
echo "##################################################"
echo "Transpose phenotype file."
awk -F'\t' '{
      for (i=1; i<=NF; i++) {
        a[NR,i] = $i;
      }
      max_nf = NF > max_nf ? NF : max_nf;
      max_nr = NR;
    }
    END {
      for (i=1; i<=max_nf; i++) {
        for (j=1; j<=max_nr; j++) {
          printf "%s%s", a[j,i], (j==max_nr ? "\n" : "\t");
        }
      }
    }' ${TMP1} > ${TMP2}

##Select phenotypes and create matrix.
echo "##################################################"
echo "Select phenotypes and create matrix."
echo ${LIST_TRAITS} | tr ',' '\n' | while read TRAIT ; do
	if [[ $(grep -P "^${TRAIT}\t" ${TMP2} | wc -l) -gt 0 ]] ; then
   	grep -P "^${TRAIT}\t" ${TMP2} >> ${TMP3}
		echo "${TRAIT}" >> ${OUTPUTFILE3}
	fi 
done

##Create phenotype matrix.
echo "##################################################"
echo "Create phenotype matrix."
cut -f2- ${TMP3} | \
awk -F'\t' '{
      for (i=1; i<=NF; i++) {
        a[NR,i] = $i;
      }
      max_nf = NF > max_nf ? NF : max_nf;
      max_nr = NR;
    }
    END {
      for (i=1; i<=max_nf; i++) {
        for (j=1; j<=max_nr; j++) {
          printf "%s%s", a[j,i], (j==max_nr ? "\n" : "\t");
        }
      }
    }' > ${OUTPUTFILE5}

##Select covariates and create matrix.
echo "##################################################"
echo "Select covariates and create matrix."
echo ${LIST_COVARIATES} | tr ',' '\n' | while read COVARIATE ; do
	if [[ $(grep -P "^${COVARIATE}\t" ${TMP2} | wc -l) -gt 0 ]] ; then
   	grep -P "^${COVARIATE}\t" ${TMP2} >> ${TMP4}
		echo "${COVARIATE}" >> ${OUTPUTFILE4}
	fi 
done

##Create covariate matrix.
echo "##################################################"
echo "Create covariate matrix."
cut -f2- ${TMP4} | \
awk -F'\t' '{
      for (i=1; i<=NF; i++) {
        a[NR,i] = $i;
      }
      max_nf = NF > max_nf ? NF : max_nf;
      max_nr = NR;
    }
    END {
      for (i=1; i<=max_nf; i++) {
        for (j=1; j<=max_nr; j++) {
          printf "%s%s", a[j,i], (j==max_nr ? "\n" : "\t");
        }
      }
    }' > ${OUTPUTFILE6}

##Create relatedness matrix (samples ordered according to the VCF file order).
##Select target samples (rows).
echo "##################################################"
echo "Create relatedness matrix (samples ordered according to the VCF file order)."
echo "Select target samples (rows)."
head -n1 ${INPUTRELATEDNESSFILE} > ${TMP5}
cat ${OUTPUTFILE2} | while read SAMPLE ; do
    grep -P "^${SAMPLE}\t" ${INPUTRELATEDNESSFILE} >> ${TMP5}
done

##Transpose relatedness file.
echo "##################################################"
echo "Transpose relatedness file."
awk -F'\t' '{
      for (i=1; i<=NF; i++) {
        a[NR,i] = $i;
      }
      max_nf = NF > max_nf ? NF : max_nf;
      max_nr = NR;
    }
    END {
      for (i=1; i<=max_nf; i++) {
        for (j=1; j<=max_nr; j++) {
          printf "%s%s", a[j,i], (j==max_nr ? "\n" : "\t");
        }
      }
    }' ${TMP5} > ${TMP6}

##Select target samples (columns).
echo "##################################################"
echo "Select target samples (columns)."
cat ${OUTPUTFILE2} | while read SAMPLE ; do
    grep -P "^${SAMPLE}\t" ${TMP6} >> ${TMP7}
done

##Transpose relatedness file.
echo "##################################################"
echo "Create relatedness matrix."
cut -f2- ${TMP7} | \
awk -F'\t' '{
      for (i=1; i<=NF; i++) {
        a[NR,i] = $i;
      }
      max_nf = NF > max_nf ? NF : max_nf;
      max_nr = NR;
    }
    END {
      for (i=1; i<=max_nf; i++) {
        for (j=1; j<=max_nr; j++) {
          printf "%s%s", a[j,i], (j==max_nr ? "\n" : "\t");
        }
      }
    }' > ${OUTPUTFILE7}

##Filter VCF file and convert it to PLINK format.
echo "##################################################"
echo "Filter VCF file and convert it to PLINK format."
CHRSETN_GAL7B=39
#plink --vcf ${OUTPUTFILE1} --maf 0.05 --geno 0.05 --allow-extra-chr --chr-set ${CHRSETN_GAL7B} --make-bed --out ${TMP_PLINKFILEPREFIX}
plink --vcf ${OUTPUTFILE1} --allow-extra-chr --chr-set ${CHRSETN_GAL7B} --make-bed --out ${TMP_PLINKFILEPREFIX}

##Run GEMMA GWAS.
echo "##################################################"
echo "Run GEMMA GWAS."
seq 1 $(cat ${OUTPUTFILE3} | wc -l) | while read PHENOTYPEi ; do
  echo "##################################################"
	echo "Running GEMMA on phenotype ${PHENOTYPEi}."
	
	##Run association analysis.
	OUTPUTFILEPREFIX=$(echo "${INPUTFILEPREFIX}.gemma_gwasoutput_phenotype${PHENOTYPEi}-job${JOBID}")
	gemma \
		-bfile ${TMP_PLINKFILEPREFIX} \
		-p ${OUTPUTFILE5} \
		-n ${PHENOTYPEi} \
		-c ${OUTPUTFILE6} \
		-k ${OUTPUTFILE7} \
		-lmm 4 \
		-outdir ${OUTPUTDIR} \
		-o ${OUTPUTFILEPREFIX}

    #Enable checks: -check
    #Disable floating point tests: -no-check 

	##Check for errors and, if no errors in output file, calculate FDR p-value.
  FILEGWAS=$(readlink -f "${OUTPUTDIR}/${OUTPUTFILEPREFIX}.assoc.txt")
	if [[ ! -s "${FILEGWAS}" ]] ; then

		echo -e "ERROR:\t${FILEGWAS}"
		echo -e "${FILEGWAS}" >> ${ERRORFILE}
		echo -e "${FILEGWAS}\tERROR" >> ${OUTPUTLIST}

	else

		echo -e "${FILEGWAS}\tCLEAR" >> ${OUTPUTLIST}
    ##Calculate FDR p-values.
    echo "Calculate FDR p-values for phenotype ${PHENOTYPEi}."
    #Select target column containing raw p-values.
    #Add original order of input p-values, sort and add column with order of ascending raw p-values.
    #Calculate BH-FDR p-values.
    #Sort p-values according to original order.
    FILEGWAS="${OUTPUTDIR}/${OUTPUTFILEPREFIX}.assoc.txt"
    NTESTS=$(cat ${FILEGWAS} | tail -n+2 | wc -l)
    cat ${FILEGWAS} | cut -f14 \
    | awk -v OFS='\t' 'NR==1 { print $0, "OriginalOrder"; next } { print $0, NR-1 }' \
    | tail -n+2 \
    | sort -g -k1,1 --temporary-directory=${TMPDIR1} \
    | awk -v OFS='\t' '{ print $0, NR }' \
    | tac | awk -v OFS='\t' -v NTESTS="${NTESTS}" '
    {
        q = ($1 * NTESTS) / $3
        if (NR == 1 || q < prev_q) prev_q = q
        if (prev_q > 1) prev_q = 1
        print $0 OFS prev_q
    }
    ' \
    | cut -f2,4 | sort -g -k1,1 --temporary-directory=${TMPDIR2} | sed '1iOriginalOrder\tPvalueFDR' > ${TMP8}
    #Save data.
    paste ${FILEGWAS} <(cut -f2 ${TMP8}) > ${TMP9}
    mv ${TMP9} ${FILEGWAS}
    #Get list of significant hits.
    TRAITNAME=$(awk -v PHENOTYPEINDEX="${PHENOTYPEi}" 'NR==PHENOTYPEINDEX' ${OUTPUTFILE3})
    OUTPUT_SIGNIFICANT="${OUTPUTDIR}/significantFDRle02-${TRAITNAME}.gemmagwas-job${JOBID}.txt"
    echo -e "Phenotype\tChromosome\tPosition\tBeta\tSE\tPvalueLRT\tPvalueFDR" > ${OUTPUT_SIGNIFICANT}
    cat ${FILEGWAS} | tail -n+2 \
    | awk -v TRAITNAME="${TRAITNAME}" '($16 <= 0.2) {print TRAITNAME"\t"$1"\t"$3"\t"$8"\t"$9"\t"$14"\t"$16}' >> ${OUTPUT_SIGNIFICANT}

	fi

done

##Remove temporary files.
echo "##################################################"
echo "Delete temporary files."
rm -rf ${TMP1} ${TMP2} ${TMP3} ${TMP4} ${TMP5} ${TMP6} ${TMP7} ${TMP8} ${TMP9} ${TMP_PLINKFILEPREFIX}* ${TMPDIR1} ${TMPDIR2}

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
Input file: ${INPUTRELATEDNESSFILE}
Input sample list: ${LIST_SAMPLES}
Input trait list: ${LIST_TRAITS}
Input covariate list: ${LIST_COVARIATES}
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