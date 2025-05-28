#!/usr/bin/env Rscript
rm(list=ls()) #Clear environment.
ARGS=commandArgs(trailingOnly=TRUE)
RUNDATE=format(Sys.time(),"%Y%m%d%H%M%S")
SCRIPTNAME<-"vcf-iHH-Rehh.R"
JOBID<-ARGS[1]
OUTPUTLOCATION<-normalizePath(ARGS[2])
sink(paste0(OUTPUTLOCATION,"/job",JOBID,".log"),type=c("output","message"))
############################################################################
################################ R SCRIPT ##################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################
#SCRIPT DESCRIPTION:

#Description:
#Run the first step of the extended haplotype homozygosity (EHH) analysis on a VCF file using Rehh.

#Input $1: Job ID.
#Input $2: Output location.
#Input $3: VCF file (.vcf.gz) split by chromosome.
#Input $4: Phased? (TRUE/FALSE).
#Output: iHH data file (.rds).

#Usage: 
#Rscript --vanilla vcf-iHH-Rehh.R <JOB ID> <OUTPUT LOCATION> <INPUT VCF FILE> <TRUE/FALSE>

############################################################################
#DIAGNOSTICS:

cat(paste0(format(Sys.time(),"%Y-%m-%d @ %H:%M:%S\n")))
cat(paste0("Submission: Rscript --vanilla ",SCRIPTNAME," ",paste(ARGS,collapse=" "),"\n"))
cat(paste0("User: ",Sys.getenv("USER"),"\n"))
cat(paste0("Home: ",Sys.getenv("HOME"),"\n"))
cat(paste0("PATH: ",Sys.getenv("PATH"),"\n"))
cat(paste0("Job ID: ",JOBID,"\n"))
sessionInfo()
cat("\n")

############################################################################
#LOAD TOOLS:

library(rehh)

############################################################################
#INPUT:

INPUTFILE<-normalizePath(ARGS[3])
INPUTFILENAME<-basename(INPUTFILE)
PHASINGSTATUS<-ARGS[4]

############################################################################
#OUTPUT:

OUTPUTFILEPREFIX<-system(paste0("echo ",INPUTFILENAME,"| sed 's/-job[0-9].*$//'"),
                         intern=TRUE)
if(is.na(PHASINGSTATUS)){PHASINGSTATUS<-FALSE} #If phasing status is not provided, assume unphased data.
if(PHASINGSTATUS){
   OUTPUTFILENAME<-paste0(OUTPUTFILEPREFIX,".phasediHH-job",JOBID,".rds")
}else{
   OUTPUTFILENAME<-paste0(OUTPUTFILEPREFIX,".unphasediHH-job",JOBID,".rds")
}
OUTPUTFILE<-paste0(OUTPUTLOCATION,"/",OUTPUTFILENAME)

############################################################################
#ACTIONS:

#Import VCF file.
cat("Converting data into an object of class haplohh (rehh::data2haplohh).\n")
HAPLODATA<-data2haplohh(hap_file=INPUTFILE,
                        polarize_vcf=FALSE, #Unpolarized data.
                        min_maf=0.05, #Filter data on a minor allele frequency or MAF.
                        vcf_reader="data.table",
                        verbose=TRUE)
#saveRDS(HAPLODATA,file=paste0(OUTPUTFILEPREFIX,".haplohh-job",JOBID,".rds"))

cat("Computing EHH based statistics over a whole chromosome (rehh::scan_hh).\n")
#Calculate iHH statistics.
HAPLODATA_iHH<-scan_hh(HAPLODATA,
                        phased=PHASINGSTATUS, #Phased data?
                        polarized=FALSE) #Unpolarized data.

#Save data.
saveRDS(HAPLODATA_iHH,file=OUTPUTFILE)

############################################################################
#SAVE CONTROL FILES:

cat(paste0("############################################################################
Date: ",RUNDATE,"
Job ID: ",JOBID,"
Script: ",SCRIPTNAME,"
Input VCF file: ",INPUTFILE,"
Phasing status: ",PHASINGSTATUS,"
Output file: ",OUTPUTFILE,"

"),file=paste0(OUTPUTLOCATION,"/README.txt"),append=TRUE)

sink()