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
##Filter out SNPs with low frequency minor genotype (SNPs that have at least of of the genotypes below teh specified frequency threshold) and index new VCF file.

##Input $1: Output location.
##Input $2: VCF file (.vcf or .vcf.gz).
##Input $3: Minimum genotype frequency (<1) or count (>1).
##Output: VCF file (.vcf.gz) and index file (.tbi).

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
##    -J vcf-filter-genotypefreq-${<INPUTFILE>##*/} \
##    vcf-filter-genotypefreq.sh <OUTPUT LOCATION> <INPUT VCF FILE> <GENOTYPE FREQUENCY THRESHOLD>

############################################################################
##SCRIPT NAME:

SCRIPTNAME=$(echo "vcf-filter-genotypefreq.sh") 

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

module load PDC/24.11
module load R/4.4.2-cpeGNU-24.11
module load bioinfo-tools
module load bcftools

############################################################################
##INPUT:

OUTPUTLOCATION=$(readlink -f $1)
INPUTFILE=$(readlink -f $2)
INPUTMINGENFREQ=$3

INPUTFILELOCATION=${INPUTFILE%/*}
INPUTFILENAME=${INPUTFILE##*/}

############################################################################
##OUTPUT:

OUTPUTFILEPREFIX=$(echo ${INPUTFILENAME} | sed 's/\.vcf.*$//' | sed 's/\.bcf.*$//' | sed 's/-job[0-9].*$//')
OUTPUTFILE1NAME=$(echo "list-SNPltMINGENFREQ.filterLFMG-job${JOBID}.txt") 
OUTPUTFILE2NAME=$(echo "list-SNPgeMINGENFREQ.filterLFMG-job${JOBID}.txt") 
OUTPUTFILE3NAME=$(echo "${OUTPUTFILEPREFIX}.filterLFMG-job${JOBID}.vcf.gz") 
OUTPUTFILE4NAME=$(echo "${OUTPUTFILE3NAME}.tbi")
OUTPUTFILE1=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE1NAME}") 
OUTPUTFILE2=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE2NAME}") 
OUTPUTFILE3=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE3NAME}") 
OUTPUTFILE4=$(echo "${OUTPUTLOCATION}/${OUTPUTFILE4NAME}") 

TMPFILE=$(echo "${OUTPUTLOCATION}/tmp-job${JOBID}.txt") 

############################################################################
##ACTIONS:

##Get matrix of translated genotypes.
{
	echo -ne "SNP\t$(bcftools query -l "${INPUTFILE}" | tr '\n' '\t' | sed 's/\t$/\n/')"
    echo ""
} > "${TMPFILE}"
bcftools query -f '%CHROM:%POS[\t%TGT]\n' "${INPUTFILE}" >> "${TMPFILE}"

##Run R to get list of SNPs with low frequency minor genotypes.
Rscript --vanilla -e '
#Read input data.
MATRIX_GENOTYPES<-read.table("'${TMPFILE}'",header=TRUE,sep="\t",stringsAsFactors=FALSE)

#Transform phased genotypes into unphased genotypes.
DATA_GENOTYPES<-data.frame(lapply(MATRIX_GENOTYPES,function(x){
x<-gsub("\\|","/",x)
x<-gsub("C/A","A/C",x)
x<-gsub("G/A","A/G",x)
x<-gsub("T/A","A/T",x)
x<-gsub("G/C","C/G",x)
x<-gsub("T/C","C/T",x)
gsub("T/G","G/T",x)}))
rownames(DATA_GENOTYPES)<-DATA_GENOTYPES[,1]

#Count samples per genotype and remove SNPs with low genotype counts (<MINGENFREQ per genotype).
GENOTYPES_COUNTS<-data.frame(SNP=DATA_GENOTYPES[,1],
                           AA=rowSums(DATA_GENOTYPES=="A/A"),
                           AC=rowSums(DATA_GENOTYPES=="A/C"),
                           AG=rowSums(DATA_GENOTYPES=="A/G"),
                           AT=rowSums(DATA_GENOTYPES=="A/T"),
                           CC=rowSums(DATA_GENOTYPES=="C/C"),
                           CG=rowSums(DATA_GENOTYPES=="C/G"),
                           CT=rowSums(DATA_GENOTYPES=="C/T"),
                           GG=rowSums(DATA_GENOTYPES=="G/G"),
                           GT=rowSums(DATA_GENOTYPES=="G/T"),
                           TT=rowSums(DATA_GENOTYPES=="T/T"),
                       	   MISSING=rowSums(DATA_GENOTYPES=="./."))
GENOTYPES_COUNTS$N_Samples<-rowSums(GENOTYPES_COUNTS[,2:12])
GENOTYPES_COUNTS$N_Genotypes<-rowSums(GENOTYPES_COUNTS[,2:12]>0)

if('${INPUTMINGENFREQ}'<1){
	GENOTYPES_COUNTS$N_Genotypes_ltMINGENFREQ<-rowSums(GENOTYPES_COUNTS[,2:12]>0 & (GENOTYPES_COUNTS[,2:12]/GENOTYPES_COUNTS$N_Samples)<MINGENFREQ)
	#SNPs that have at least one genotype with frequency <MINGENFREQ.
	SNPltMINGENFREQ<-GENOTYPES_COUNTS$SNP[which(rowSums(GENOTYPES_COUNTS[,2:12]>0 & (GENOTYPES_COUNTS[,2:12]/GENOTYPES_COUNTS$N_Samples)<MINGENFREQ)>0)]
	SNPgeMINGENFREQ<-GENOTYPES_COUNTS$SNP[!(GENOTYPES_COUNTS$SNP %in% SNPltMINGENFREQ)]
}else{
	GENOTYPES_COUNTS$N_Genotypes_ltMINGENFREQ<-rowSums(GENOTYPES_COUNTS[,2:12]>0 & GENOTYPES_COUNTS[,2:12]<'${INPUTMINGENFREQ}')
	#SNPs that have at least one genotype with frequency <MINGENFREQ.
	SNPltMINGENFREQ<-GENOTYPES_COUNTS$SNP[which(rowSums(GENOTYPES_COUNTS[,2:12]>0 & GENOTYPES_COUNTS[,2:12]<'${INPUTMINGENFREQ}')>0)]
	SNPgeMINGENFREQ<-GENOTYPES_COUNTS$SNP[!(GENOTYPES_COUNTS$SNP %in% SNPltMINGENFREQ)]
}

write.table(matrix(SNPltMINGENFREQ),"'${OUTPUTFILE1}'",sep="\n",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(matrix(SNPgeMINGENFREQ),"'${OUTPUTFILE2}'",sep="\n",row.names=FALSE,col.names=FALSE,quote=FALSE)
'

#Filter VCF file.
bcftools view --regions $(cat ${OUTPUTFILE2} | paste -sd,) ${INPUTFILE} | bcftools sort --temp-dir "${OUTPUTLOCATION}/" --write-index=tbi --output-type z --output ${OUTPUTFILE3}

rm -f ${TMPFILE}

############################################################################
##SAVE CONTROL FILES:

if [[ -s "${OUTPUTFILE4}" ]]; then
	
	##README LOG ENTRY:
echo "############################################################################
Date: ${RUNDATE}
Job ID: ${JOBID}
Script: ${PATHTOMYSUBMITTEDSCRIPTS}/job${JOBID}-date${RUNDATE}-${SCRIPTNAME}
Input file: ${INPUTFILE}
Input minimum genotype frequency: ${INPUTGENFREQ}
Output file: ${OUTPUTFILE1}
Output file: ${OUTPUTFILE2}
Output file: ${OUTPUTFILE3}
Output file: ${OUTPUTFILE4}
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