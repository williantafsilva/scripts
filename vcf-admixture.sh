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
##Align FASTA sequences using ClustalOmega.

##Input $1: Output location.
##Input $2: Indexed VCF file (.vcf.gz).
##Input $3: Comma-separated list of K values (possible numbers of ancestral populations; use paste -sd,).
##Input $4: Comma-separated list of samples to be included in the analysis. If not provided, all samples will be used.
##Output: Admixture files (*.Q, *.P, *.Q_se).

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
##    -J vcf-admixture-${<INPUTFILE>##*/} \
##    vcf-admixture.sh <OUTPUT LOCATION> <INPUT VCF FILE> <K> <SAMPLES>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-admixture.sh")

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
module load ADMIXTURE/1.3.0

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)
INPUTLIST_K=$3
INPUTLIST_SAMPLES=$4

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

INPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/-job[0-9].*$//')

TMPFILEPREFIX="${OUTPUTLOCATION}/tmp-job${JOBID}"

############################################################################
##ACTIONS:

#Count chromosomes.
CHRSET=$(bcftools index -s ${INPUTFILE} | cut -f 1 | wc -l)

#Get list of samples.
if [[ -z "${INPUTLIST_SAMPLES}" ]]; then
	SAMPLELIST=$(bcftools query -l ${INPUTFILE} | paste -sd,)
else
	SAMPLELIST=${INPUTLIST_SAMPLES}
fi

#Select target samples and filter VCF file, keeping only biallelic SNPs.
bcftools view \
--samples ${SAMPLELIST} \
--types snps \
--min-alleles 2 \
--max-alleles 2 \
--min-ac 1 \
--write-index=tbi \
--output-type z \
--output "${TMPFILEPREFIX}.vcf.gz" \
${INPUTFILE}

#Convert VCF to PLINK format.
plink --vcf "${TMPFILEPREFIX}.vcf.gz" --chr-set ${CHRSET} --allow-extra-chr --make-bed --out ${TMPFILEPREFIX}

#Perform LD pruning.
#Prune according to a correlation threshold and store the pruned dataset in prunedData.bed
#(remove each SNP that has an R2 value greater than 0.1 with any other SNP within a 50-SNP 
#sliding window (advanced by 10 SNPs each time).
plink --bfile ${TMPFILEPREFIX} --chr-set ${CHRSET} --allow-extra-chr --indep-pairwise 50 10 0.1 --out "${TMPFILEPREFIX}_pruned"
plink --bfile ${TMPFILEPREFIX} --chr-set ${CHRSET} --allow-extra-chr --extract "${TMPFILEPREFIX}_pruned.prune.in" --make-bed --out "${TMPFILEPREFIX}_filtered"

#Run ADMIXTURE.
for K in $(tr ',' '\n' <<< "${INPUTLIST_K}"); do 

	#Create output file prefix for each K.
	OUTPUTFILEPREFIX="${OUTPUTLOCATION}/${INPUTFILEPREFIX}.admixtureK${K}-job${JOBID}"

	#Run ADMIXTURE.
	admixture -B100 --cv=5 "${TMPFILEPREFIX}_filtered.bed" ${K} -j10 | tee "${OUTPUTFILEPREFIX}.out"

	#Ancestral population names.
	POPNAMES=""
	for ((i=1; i<=K; i++)); do
    	POPNAMES+="AncestralPopulation${i}"$'\t'
	done
	#Remove the trailing tab.
	POPNAMES=${POPNAMES%$'\t'}

	#Add row names and header.
	#Allele frequencies of the inferred ancestral populations.
	paste <(bcftools view -H "${TMPFILEPREFIX}.vcf.gz" | cut -f1 | sed '1iChromosome') \
	<(bcftools view -H "${TMPFILEPREFIX}.vcf.gz" | cut -f2 | sed '1iPosition') \
	<(cat "${TMPFILEPREFIX}_filtered.${K}.P" | tr ' ' '\t' | sed "1i${POPNAMES}") > "${OUTPUTFILEPREFIX}.P"
	#Point-estimate admixture proportions (ancestry fractions).
	paste <(bcftools query -l "${TMPFILEPREFIX}.vcf.gz" | sed '1iSample') \
	<(cat "${TMPFILEPREFIX}_filtered.${K}.Q" | tr ' ' '\t' | sed "1i${POPNAMES}") > "${OUTPUTFILEPREFIX}.Q"
	#Bootstrap standard errors (se) for each of the ancestry proportion estimates.
	paste <(bcftools query -l "${TMPFILEPREFIX}.vcf.gz" | sed '1iSample') \
	<(cat "${TMPFILEPREFIX}_filtered.${K}.Q_se" | tr ' ' '\t' | sed "1i${POPNAMES}") > "${OUTPUTFILEPREFIX}.Q_se"
	#Estimated statistical bias of the individual ancestry proportions (the Q values) for each of the sampled individuals.
	paste <(bcftools query -l "${TMPFILEPREFIX}.vcf.gz" | sed '1iSample') \
	<(cat "${TMPFILEPREFIX}_filtered.${K}.Q_bias" | tr ' ' '\t' | sed "1i${POPNAMES}") > "${OUTPUTFILEPREFIX}.Q_bias"

done

#Delete temporary files.
rm -f $(find ${TMPFILEPREFIX}*)

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTLOCATION}/${INPUTFILEPREFIX}.admixtureK${K}-job${JOBID}.Q_se" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input VCF file: ${INPUTFILE}
Input list of K values: ${INPUTLIST_K}
Input list of samples: ${INPUTLIST_SAMPLES}
$(readlink -f $(find ${OUTPUTLOCATION}/*-job${JOBID}* -maxdepth 1) | sed 's/^/Output file: /g')
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