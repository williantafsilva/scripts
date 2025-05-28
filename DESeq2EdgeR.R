#!/usr/bin/env Rscript
rm(list=ls()) #Clear environment.
ARGS=commandArgs(trailingOnly=TRUE) 
RUNDATE=format(Sys.time(),"%Y%m%d%H%M%S")
SCRIPTNAME<-"DESeq2EdgeR.R"
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
#Run differential expression analyses using DESeq2 and EdgeR.

#Input $1: Job ID.
#Input $2: Output location.
#Input $3: Expression matrix file (tab-separated).
#Input $4: Metadata/groups file (tab-separated).
#Input $5: Formula (e.g., ~Sex+Age).
#Output: DESeq2 and EdgeR results (*.txt) and volcano plots (*.pdf).

#Usage: 
#Rscript --vanilla DESeq2EdgeR.R <JOB ID> <OUTPUT LOCATION> <EXPRESSION FILE> <METADATA FILE> <FORMULA>

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
library(ggfortify)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(tibble)

############################################################################
#INPUT:

#Count data file.
FILE_COUNTDATA<-normalizePath(ARGS[3])
DATA_COUNTS<-read.table(FILE_COUNTDATA,header=TRUE,row.names=1,sep="\t",stringsAsFactors=FALSE)

#Metadata/groups file.
FILE_METADATA<-normalizePath(ARGS[4])
DATA_METADATA<-read.table(FILE_METADATA,header=TRUE,row.names=1,sep="\t",stringsAsFactors=TRUE)

#Formula.
DESIGNFORMULA=as.formula(ARGS[5])

############################################################################
#OUTPUT:

OUTPUTFILE1NAME<-paste0("DESeq2result-job",JOBID,".txt")
OUTPUTFILE2NAME<-paste0("EdgeRresult-job",JOBID,".txt")
OUTPUTFILE3NAME<-paste0("DESeq2plot-job",JOBID,".pdf")
OUTPUTFILE4NAME<-paste0("EdgeRplot-job",JOBID,".pdf")
OUTPUTFILE1<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE1NAME)
OUTPUTFILE2<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE2NAME)
OUTPUTFILE3<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE3NAME)
OUTPUTFILE4<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE4NAME)

############################################################################
#ACTIONS:

#Ensure that sample names are the same in counts
if(all(colnames(DATA_COUNTS) %in% colnames(DATA_METADATA))){
  
  tDATA_METADATA<-t(DATA_METADATA)
  
  ######################################
  #Run DESeq2 analysis.
  DESeq2_DDS<-DESeqDataSetFromMatrix(countData=round(DATA_COUNTS),colData=tDATA_METADATA,design=DESIGNFORMULA)
  DESeq2_DDS<-DESeq(DESeq2_DDS)
  DESeq2_RESULTS<-results(DESeq2_DDS)
  DESeq2_RESULTS<-rownames_to_column(as.data.frame(DESeq2_RESULTS),var="FeatureID")
  
  #Plot DESeq2 results.
  p.DESeq2<-RESULTS_DESeq2 %>% 
    filter(!is.na(padj)) %>%
    mutate(Color=ifelse(padj<0.01,"blue","black")) %>%
    mutate(Color=ifelse(padj<0.01 & abs(log2FoldChange)>2,"red",Color)) %>%
    ggplot(aes(x=log2FoldChange,y=-log10(padj),color=Color))+
    geom_point(alpha=1,
               shape=20,
               size=1)+
    geom_blank(aes(x=-log2FoldChange))+
    geom_hline(yintercept=0,color="black",linetype="dashed")+
    labs(x="Log2FoldChange",y="-Log10(Pvalue)")+
    ggtitle("DESeq2")+
    myggplottheme_nolegend
  
  ######################################
  #Rund edgeR analysis.
  EDGER_DGE<-DGEList(counts=round(DATA_COUNTS))
  EDGER_DGE<-calcNormFactors(EDGER_DGE)
  EDGER_DESIGN<-model.matrix(DESIGNFORMULA,data=as.data.frame(tMETADATA))
  EDGER_DGE<-estimateDisp(EDGER_DGE,EDGER_DESIGN)
  EDGER_FIT<-glmFit(EDGER_DGE,EDGER_DESIGN)
  EDGER_LRT<-glmLRT(EDGER_FIT)
  EDGER_RESULTS<-topTags(EDGER_LRT,n=Inf)$table
  EDGER_RESULTS<-rownames_to_column(EDGER_RESULTS,var="FeatureID")
  
  #Plot edgeR results.
  p.EdgeR<-EDGER_RESULTS %>% 
    filter(!is.na(FDR)) %>%
    mutate(Color=ifelse(FDR<0.01,"blue","black")) %>%
    mutate(Color=ifelse(FDR<0.01 & abs(logFC)>2,"red",Color)) %>%
    ggplot(aes(x=logFC,y=-log10(FDR),color=Color))+
    geom_point(alpha=1,
               shape=20,
               size=1)+
    geom_blank(aes(x=-logFC))+
    geom_hline(yintercept=0,color="black",linetype="dashed")+
    labs(x="Log2FoldChange",y="-Log10(Pvalue)")+
    ggtitle("EdgeR")+
    myggplottheme_nolegend

  ######################################
  #Save results.
  write.table(DESeq2_RESULTS,file=OUTPUTFILE1,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE,append=FALSE)
  write.table(EDGER_RESULTS,file=OUTPUTFILE2,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE,append=FALSE)
  ggsave(OUTPUTFILE3,plot=p.DESeq2,width=15,height=15,units="cm")
  ggsave(OUTPUTFILE4,plot=p.EdgeR,width=15,height=15,units="cm")
  
}else{
  stop("Sample names in counts data do not match sample names in metadata.")
}

############################################################################
#SAVE CONTROL FILES:

cat(paste0("############################################################################
Date: ",RUNDATE,"
Job ID: ",JOBID,"
Script: ",SCRIPTNAME,"
Input count file: ",FILE_COUNTDATA,"
Input metadata file: ",FILE_METADATA,"
Input formula: ",ARGS[5],"
Output file: ",OUTPUTFILE1,"
Output file: ",OUTPUTFILE2,"
Output file: ",OUTPUTFILE3,"
Output file: ",OUTPUTFILE4,"

"),file=paste0(OUTPUTLOCATION,"/README.txt"),append=TRUE)

sink()