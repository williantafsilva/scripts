#!/usr/bin/env Rscript
rm(list=ls()) #Clear environment.
ARGS=commandArgs(trailingOnly=TRUE)
RUNDATE=format(Sys.time(),"%Y%m%d%H%M%S")
SCRIPTNAME<-"iHH-iHS-Rehh.R"
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
#Run the second step of the extended haplotype homozygosity (EHH) analysis on iHH data calculated in the first step of the analysis.

#Input $1: Job ID.
#Input $2: Output location.
#Input $3: Directory containing iHH data files (.rds) split by chromosome.
#Output: iHS data file (.rds).

#Usage: 
#Rscript --vanilla iHH-iHS-Rehh.R <JOB ID> <OUTPUT LOCATION> <INPUT DIRECTORY>

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

library(R.utils)
library(tidyverse)
library(rehh)

############################################################################
#INPUT:

INPUTDIR<-normalizePath(ARGS[3])
INPUTDIRNAME<-basename(INPUTDIR)

############################################################################
#OUTPUT:

OUTPUTFILEPREFIX<-system(paste0("echo ",INPUTDIRNAME,"| sed 's/-job[0-9].*$//'"),
                         intern=TRUE)
OUTPUTFILE1NAME<-paste0(OUTPUTFILEPREFIX,".RehhiHS-job",JOBID,".txt")
OUTPUTFILE2NAME<-paste0(OUTPUTFILEPREFIX,".RehhiHSfreq-job",JOBID,".txt")
OUTPUTFILE3NAME<-paste0(OUTPUTFILEPREFIX,".RehhiHSregions-job",JOBID,".txt")
OUTPUTFILE1<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE1NAME)
OUTPUTFILE2<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE2NAME)
OUTPUTFILE3<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE3NAME)

############################################################################
#ACTIONS:

#Set working directory.
setwd(INPUTDIR)

#List of files in the input directory.
FILES<-list.files(INPUTDIR) 
FILES<-FILES[grep("iHH-job",FILES)]
FILES<-FILES[grep(".rds",FILES)]

cat(paste0("Processing file: ",FILES[1],".\n"))
iHHFILE<-normalizePath(FILES[1])
iHH_DATA<-readRDS(iHHFILE)
for(C in 2:length(FILES)){ #for each chromosome.
  cat(paste0("Processing file: ",FILES[C],".\n"))
  iHHFILE<-normalizePath(FILES[C])
  iHH_DATA_TMP<-readRDS(iHHFILE)
  iHH_DATA<-rbind(iHH_DATA,iHH_DATA_TMP)
}

#Calculate genome-wide iHS values
iHS_DATA<-ihh2ihs(iHH_DATA,
                  freqbin=0.025,
                  p.adjust.method="BH")

#Calculate candidate regions.
CANDIDATEREGIONS<-calc_candidate_regions(iHS_DATA,
                                         threshold=2,
                                         pval=TRUE,
                                         window_size=20000,
                                         overlap=2000,
                                         min_n_extr_mrk=2)

#Change column names.
colnames(iHS_DATA$ihs)<-c("Chromosome","Position","iHS","LogPvalue")

#Remove NAs and sort by chromosome and position.
iHS_DATA$ihs<-iHS_DATA$ihs %>%
  filter(!is.na(LogPvalue)) %>%
  arrange(Chromosome,Position)

#Save data.
write.table(iHS_DATA$ihs,OUTPUTFILE1,sep="\t",row.names=FALSE,quote=FALSE)
write.table(iHS_DATA$frequency.class,OUTPUTFILE2,sep="\t",row.names=FALSE,quote=FALSE)
write.table(CANDIDATEREGIONS,OUTPUTFILE3,sep="\t",row.names=FALSE,quote=FALSE)

############################################################################
#SAVE CONTROL FILES:

cat(paste0("############################################################################
Date: ",RUNDATE,"
Job ID: ",JOBID,"
Script: ",SCRIPTNAME,"
Input directory: ",INPUTDIR,"
Output file: ",OUTPUTFILE1,"
Output file: ",OUTPUTFILE2,"
Output file: ",OUTPUTFILE3,"

"),file=paste0(OUTPUTLOCATION,"/README.txt"),append=TRUE)

sink()