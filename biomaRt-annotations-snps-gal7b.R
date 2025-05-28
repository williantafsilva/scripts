#!/usr/bin/env Rscript
rm(list=ls()) #Clear environment.
ARGS=commandArgs(trailingOnly=TRUE)
RUNDATE=format(Sys.time(),"%Y%m%d%H%M%S")
SCRIPTNAME<-"biomaRt-annotations-snps-gal7b.R"
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
#Get SNP annotations for a genomic region, given the chromosome, the start position 
#and the end position of the region of interest.

#Input $1: Job ID.
#Input $2: Output location.
#Input $3: Output file tag.
#Input $4: Chromosome (or comma-separated list of chromosome names; use | paste -sd,).
#Input $5: Start position (or comma-separated list of start positions; use | paste -sd,).
#Input $6: End position (or comma-separated list of end positions; use | paste -sd,).
#Output: File with annotations (*.txt).

#Usage: 
#Rscript --vanilla biomaRt-annotations-snps-gal7b.R <JOB ID> <OUTPUT LOCATION> <FILE TAG> <CHROMOSOME> <START POSITION> <END POSITIONS> 

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

library(tidyverse)
library(biomaRt)

############################################################################
#INPUT:

OUTPUTFILETAG<-ARGS[3]
CHRLIST<-ARGS[4]
STARTPOSLIST<-ARGS[5]
ENDPOSLIST<-ARGS[6]

############################################################################
#OUTPUT:

OUTPUTFILENAME<-paste0("biomartgal7b.annotsnps",OUTPUTFILETAG,"-job",JOBID,".txt")
OUTPUTFILE<-paste0(OUTPUTLOCATION,"/",OUTPUTFILENAME)

############################################################################
#ACTIONS:

#Input data.
INPUTDATA<-data.frame(Chromosome=strsplit(CHRLIST,",")[[1]],
                      Start=strsplit(STARTPOSLIST,",")[[1]],
                      End=strsplit(ENDPOSLIST,",")[[1]])
INPUTDATA$Chromosome<-gsub("Chr","",INPUTDATA$Chromosome) #Remove "Chr" from chromosome names.
INPUTDATA$Chromosome<-gsub("chr","",INPUTDATA$Chromosome) #Remove "chr" from chromosome names.
INPUTDATA

NREGIONS<-nrow(INPUTDATA)
cat(paste0("Number of input regions: ",NREGIONS,"\n"))

#Clear cache.
#biomartCacheClear()

#Load Mart object with desired dataset.
ensembl_snps_gal7b<-useEnsembl(biomart="snps",dataset="ggallus_snp") #SNPs.

#Query database.
cat(paste0("Reading region 1: ",INPUTDATA$Chromosome[1],":",INPUTDATA$Start[1],"-",INPUTDATA$End[1],"\n"))
DATA_SNPS<-getBM(attributes=c("chr_name","chrom_start","refsnp_id",
                              "allele","ensembl_gene_stable_id","ensembl_type"),
                 filters=c("chr_name","start","end"),
                 values=list("chr_name"=INPUTDATA$Chromosome[1],
                             "start"=INPUTDATA$Start[1],
                             "end"=INPUTDATA$End[1]),
                 mart=ensembl_snps_gal7b)

if(nrow(DATA_SNPS)==0){
  DATA_SNPS<-data.frame(chr_name=NA,
                        chrom_start=NA,
                        refsnp_id=NA,
                        allele=NA,
                        ensembl_gene_stable_id=NA,
                        ensembl_type=NA,
                        InputRegion=paste0(INPUTDATA$Chromosome[1],":",INPUTDATA$Start[1],"-",INPUTDATA$End[1]))
}else{
  DATA_SNPS[DATA_SNPS==""]<-NA
  DATA_SNPS$InputRegion<-paste0(INPUTDATA$Chromosome[1],":",INPUTDATA$Start[1],"-",INPUTDATA$End[1])
}

if(nrow(INPUTDATA)>1){
  for(i in 2:nrow(INPUTDATA)){
    cat(paste0("Reading region ",i,"/",NREGIONS,": ",INPUTDATA$Chromosome[i],":",INPUTDATA$Start[i],"-",INPUTDATA$End[i],"\n"))
    TMP<-getBM(attributes=c("chr_name","chrom_start","refsnp_id",
                              "allele","ensembl_gene_stable_id","ensembl_type"),
                 filters=c("chr_name","start","end"),
                 values=list("chr_name"=INPUTDATA$Chromosome[i],
                             "start"=INPUTDATA$Start[i],
                             "end"=INPUTDATA$End[i]),
                 mart=ensembl_snps_gal7b)
    
    if(nrow(TMP)>0){
      TMP[TMP==""]<-NA
      TMP$InputRegion<-paste0(INPUTDATA$Chromosome[i],":",INPUTDATA$Start[i],"-",INPUTDATA$End[i])
      DATA_SNPS<-rbind(DATA_SNPS,TMP)
    }
    cat(paste0("Number of hits: ",nrow(TMP),"\n"))
  }
}
DATA_SNPS<-DATA_SNPS[!is.na(DATA_SNPS$chr_name),]

cat(paste0("Total number of hits: ",nrow(DATA_GENES),"\n"))

#Save formatted data.
cat(paste0("Saving data."))
write.table(DATA_SNPS,file=OUTPUTFILE,sep="\t",row.names=FALSE,col.names=TRUE)

############################################################################
#SAVE CONTROL FILES:

cat(paste0("############################################################################
Date: ",RUNDATE,"
Job ID: ",JOBID,"
Script: ",SCRIPTNAME,"
Output file tag: ",OUTPUTFILETAG,"
Chromosome: ",CHRLIST,"
Start position: ",STARTPOSLIST,"
End position: ",ENDPOSLIST,"
Output file: ",OUTPUTFILE,"

"),file=paste0(OUTPUTLOCATION,"/README.txt"),append=TRUE)

sink()