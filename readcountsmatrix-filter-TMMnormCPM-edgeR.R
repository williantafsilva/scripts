#!/usr/bin/env Rscript
rm(list=ls()) #Clear environment.
ARGS=commandArgs(trailingOnly=TRUE)
RUNDATE=format(Sys.time(),"%Y%m%d%H%M%S")
SCRIPTNAME<-"readcountsmatrix-filter-TMMnormCPM-edgeR.R"
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
#Input $3: Read counts matrix (feature/gene/transcript ID per row and sample per column).
#Input $4: Plot data? (TRUE/FALSE).
#Output: TMM-normalized counts and TMM-normalized CPM files (.txt).

#Usage: 
#Rscript --vanilla readcountsmatrix-filter-TMMnorm-edgeR.R <JOB ID> <OUTPUT LOCATION> <INPUT FILE> <TRUE/FALSE>

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
library(ggplot2)
library(edgeR)
#library(DESeq2)
library(tibble)
library(ggpubr)

myggplottheme_blank<-theme(title=element_text(size=10,face="bold"),
                           axis.title=element_text(size=10,face="bold"),
                           axis.text=element_text(size=10),
                           axis.text.x=element_text(angle=60,size=8,vjust=0.5),
                           #legend.position="none",
                           legend.title=element_text(size=10,face="bold"),
                           legend.text=element_text(size=10),
                           legend.key=element_blank(),
                           panel.grid=element_line(colour="gray90"),
                           panel.grid.major.x=element_blank(),
                           panel.grid.minor.x=element_blank(),
                           panel.background=element_rect(fill="white",colour="black"),
                           panel.grid.major=element_blank(),
                           panel.grid.minor=element_blank(),
                           strip.background=element_rect(colour="black",
                                                         fill="white"))

############################################################################
#INPUT:

INPUTFILE<-normalizePath(ARGS[3])
if(is.na(ARGS[4])){PLOTDATA<-FALSE}else{PLOTDATA<-ARGS[4]}

INPUTFILENAME<-basename(INPUTFILE)

############################################################################
#OUTPUT:

OUTPUTFILEPREFIX<-system(paste0("echo ",INPUTFILENAME," | sed 's/.txt$//' | sed 's/-job[0-9].*$//'"),
                         intern=TRUE)
OUTPUTFILE1NAME<-paste0(OUTPUTFILEPREFIX,".TMMnormCPMraw-job",JOBID,".txt")
OUTPUTFILE2NAME<-paste0(OUTPUTFILEPREFIX,".TMMnormcountsraw-job",JOBID,".txt")
OUTPUTFILE3NAME<-paste0(OUTPUTFILEPREFIX,".TMMnormCPMfiltered-job",JOBID,".txt")
OUTPUTFILE4NAME<-paste0(OUTPUTFILEPREFIX,".TMMnormcountsfiltered-job",JOBID,".txt")
OUTPUTFILE5NAME<-paste0(OUTPUTFILEPREFIX,".TMMnormplotsraw-job",JOBID,".pdf")
OUTPUTFILE6NAME<-paste0(OUTPUTFILEPREFIX,".TMMnormplotsfiltered-job",JOBID,".pdf")
OUTPUTFILE1<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE1NAME)
OUTPUTFILE2<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE2NAME)
OUTPUTFILE3<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE3NAME)
OUTPUTFILE4<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE4NAME)
OUTPUTFILE5<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE5NAME)
OUTPUTFILE6<-paste0(OUTPUTLOCATION,"/",OUTPUTFILE6NAME)

############################################################################
#ACTIONS:

#Read input read counts matrix.
COUNTS<-read.table(INPUTFILE,header=TRUE,row.names=1,sep="\t",check.names=FALSE)

######################################
#Plots (pre-filtering).
if(PLOTDATA==TRUE){
  BOXDATA<-stack(as.data.frame(log2(as.matrix(COUNTS)+1)))
  p1<-ggplot(BOXDATA,aes(x=ind,y=values))+ 
    geom_boxplot(outlier.shape=16,outlier.size=0.3,
                 outlier.colour="red",outlier.fill="red")+
    xlab("Sample")+
    ylab(expression('Log'[2]~'Read Counts'))+
    ggtitle("Raw data")+
    myggplottheme_blank+
    theme(axis.text.x=element_blank())
  #p1

  HISTDATA<-data.frame(Count=c(log2(as.matrix(COUNTS)+1)))
  p2<-ggplot(HISTDATA,aes(x=Count))+
    geom_histogram(aes(y=after_stat(density)),fill="cyan",color="#e9ecef",alpha=0.9)+
    geom_density()+
    geom_vline(aes(xintercept=mean(Count)),
               color="black",linetype="dashed",linewidth=0.5)+
    xlab(expression('Log'[2]~'Read Counts'))+
    ylab("Density")+
    #ggtitle("Raw data")+
    myggplottheme_blank
  #p2

  BARGDATA<-data.frame(ColSum=colSums(COUNTS>0))
  p3<-ggplot(data=BARGDATA,aes(x=row.names(BARGDATA),y=ColSum))+
    geom_bar(stat="identity",fill="purple")+
    geom_hline(aes(yintercept=median(colSums(COUNTS>0))),
               color="black",linetype="dashed",linewidth=0.5)+
    xlab("Sample")+
    ylab("Number of detected features")+
    #ggtitle("Raw data")+
    myggplottheme_blank+
    theme(axis.text.x=element_blank())
  #p3

  BARSDATA<-data.frame(RowSum=rowSums(COUNTS>0))
  p4<-ggplot(data=BARSDATA,aes(x=row.names(BARSDATA),y=RowSum))+
    geom_bar(stat="identity",fill="orange")+
    geom_hline(aes(yintercept=median(rowSums(COUNTS>0))),
               color="red",linetype="dashed",linewidth=0.5)+
    xlab("Feature")+
    ylab("Number of samples")+
    #ggtitle("Raw data")+
    myggplottheme_blank+
    theme(axis.text.x=element_blank())
  #p4

  p5<-ggplot(BARSDATA,aes(x=RowSum))+
    geom_histogram(aes(y=after_stat(density)),fill="cyan",color="#e9ecef",alpha=0.9)+
    geom_density()+
    geom_vline(aes(xintercept=mean(RowSum)),
               color="black",linetype="dashed",linewidth=0.5)+
    xlab("Number of samples")+
    ylab("Density")+
    #ggtitle("Raw data")+
    myggplottheme_blank
  #p5

  p.raw.all<-ggarrange(p1,p2,p3,p4,p5,
                       nrow=5) 
  #p.raw.all
}

######################################
#Pre-filtering normalization.

#Remove samples with all-zero counts.
COUNTS<-COUNTS[,colSums(COUNTS)!=0]

#Create DGEList object.
RAW_DGEOBJECT<-DGEList(counts=COUNTS)

#Calculate normalization factors using TMM.
RAW_DGEOBJECT<-calcNormFactors(RAW_DGEOBJECT,method="TMM")

#Compute normalized counts (CPM with normalization factors).
RAW_NORMCPMCOUNTS<-cpm(RAW_DGEOBJECT,normalized.lib.sizes=TRUE)
RAW_NORMCPMCOUNTS<-tibble::rownames_to_column(as.data.frame(RAW_NORMCPMCOUNTS),"ID")

#Apply normalization factors to raw counts.
RAW_NORMFACTORS<-RAW_DGEOBJECT$samples$norm.factors
RAW_LIBSIZES<-RAW_DGEOBJECT$samples$lib.size
RAW_EFFECTIVELIBSIZES<-RAW_LIBSIZES*RAW_NORMFACTORS

#Scale counts.
RAW_NORMCOUNTS<-t((t(RAW_DGEOBJECT$counts)/RAW_EFFECTIVELIBSIZES)*mean(RAW_EFFECTIVELIBSIZES))
RAW_NORMCOUNTS<-tibble::rownames_to_column(as.data.frame(RAW_NORMCOUNTS),"ID")

######################################
#Filter out samples and features (e.g., genes).
NFEATURES<-max(colSums(COUNTS>0))
NSAMPLES<-max(rowSums(COUNTS>0))
MINREADS<-10 #Minimum accepted number of reads per feature.
MINFEATURES<-floor(0.10*NFEATURES) #Minimum fraction of detected features per sample.
MINSAMPLES<-floor(0.95*NSAMPLES) #Minimum number of samples per detected feature.

#Keep only features with at least MINREADS reads in at least MINSAMPLES samples.
KEEPFEATURES<-rowSums(COUNTS>MINREADS)>=MINSAMPLES
COUNTS<-COUNTS[KEEPFEATURES,]

#Keep only samples with at least MINFEATURES detected features.
KEEPSAMPLES<-colSums(COUNTS>MINREADS)>=MINFEATURES
COUNTS<-COUNTS[,KEEPSAMPLES]

######################################
if(nrow(COUNTS)>0 & ncol(COUNTS)>0){ #Only if filtered data is not empty.
  #Plots (post-filtering).
  if(PLOTDATA==TRUE){
    BOXDATA<-stack(as.data.frame(log2(as.matrix(COUNTS)+1)))
    p6<-ggplot(BOXDATA,aes(x=ind,y=values))+ 
      geom_boxplot(outlier.shape=16,outlier.size=0.3,
                   outlier.colour="red",outlier.fill="red")+
      xlab("Sample")+
      ylab(expression('Log'[2]~'Read Counts'))+
      ggtitle("Filtered data")+
      myggplottheme_blank+
      theme(axis.text.x=element_blank())
    #p6

    HISTDATA<-data.frame(Count=c(log2(as.matrix(COUNTS)+1)))
    p7<-ggplot(HISTDATA,aes(x=Count))+
      geom_histogram(aes(y=after_stat(density)),fill="cyan",color="#e9ecef",alpha=0.9)+
      geom_density()+
      geom_vline(aes(xintercept=mean(Count)),
                 color="black",linetype="dashed",linewidth=0.5)+
      xlab(expression('Log'[2]~'Read Counts'))+
      ylab("Density")+
      #ggtitle("Raw data")+
      myggplottheme_blank
    #p7

    BARGDATA<-data.frame(ColSum=colSums(COUNTS>0))
    p8<-ggplot(data=BARGDATA,aes(x=row.names(BARGDATA),y=ColSum))+
      geom_bar(stat="identity",fill="purple")+
      geom_hline(aes(yintercept=median(colSums(COUNTS>0))),
                 color="black",linetype="dashed",linewidth=0.5)+
      xlab("Sample")+
      ylab("Number of detected features")+
      #ggtitle("Raw data")+
      myggplottheme_blank+
      theme(axis.text.x=element_blank())
    #p8

    BARSDATA<-data.frame(RowSum=rowSums(COUNTS>0))
    p9<-ggplot(data=BARSDATA,aes(x=row.names(BARSDATA),y=RowSum))+
      geom_bar(stat="identity",fill="orange")+
      geom_hline(aes(yintercept=median(rowSums(COUNTS>0))),
                 color="red",linetype="dashed",linewidth=0.5)+
      xlab("Feature")+
      ylab("Number of samples")+
      #ggtitle("Raw data")+
      myggplottheme_blank+
      theme(axis.text.x=element_blank())
    #p9

    p10<-ggplot(BARSDATA,aes(x=RowSum))+
      geom_histogram(aes(y=after_stat(density)),fill="cyan",color="#e9ecef",alpha=0.9)+
      geom_density()+
      geom_vline(aes(xintercept=mean(RowSum)),
                 color="black",linetype="dashed",linewidth=0.5)+
      xlab("Number of samples")+
      ylab("Density")+
      #ggtitle("Raw data")+
      myggplottheme_blank
    #p10

    p.filtered.all<-ggarrange(p6,p7,p8,p9,p10,
                              nrow=5) 
    #p.filtered.all
  }

  ######################################
  #Normalization.
  #Create DGEList object.
  DGEOBJECT<-DGEList(counts=COUNTS)

  #Calculate normalization factors using TMM.
  DGEOBJECT<-calcNormFactors(DGEOBJECT,method="TMM")

  #Compute normalized counts (CPM with normalization factors).
  NORMCPMCOUNTS<-cpm(DGEOBJECT,normalized.lib.sizes=TRUE)
  NORMCPMCOUNTS<-tibble::rownames_to_column(as.data.frame(NORMCPMCOUNTS),"ID")

  #Apply normalization factors to raw counts.
  NORMFACTORS<-DGEOBJECT$samples$norm.factors
  LIBSIZES<-DGEOBJECT$samples$lib.size
  EFFECTIVELIBSIZES<-LIBSIZES*NORMFACTORS

  #Scale counts.
  NORMCOUNTS<-t((t(DGEOBJECT$counts)/EFFECTIVELIBSIZES)*mean(EFFECTIVELIBSIZES))
  NORMCOUNTS<-tibble::rownames_to_column(as.data.frame(NORMCOUNTS),"ID")
}

######################################
#Save data.
write.table(RAW_NORMCPMCOUNTS,OUTPUTFILE1,sep="\t",row.names=FALSE,quote=FALSE)
write.table(RAW_NORMCOUNTS,OUTPUTFILE2,sep="\t",row.names=FALSE,quote=FALSE)

#Save plots.
if(PLOTDATA==TRUE){
  ggsave(OUTPUTFILE5,p.raw.all,width=10,height=20)
}

if(nrow(COUNTS)>0 & ncol(COUNTS)>0){
  #Save data.
  write.table(NORMCPMCOUNTS,OUTPUTFILE3,sep="\t",row.names=FALSE,quote=FALSE)
  write.table(NORMCOUNTS,OUTPUTFILE4,sep="\t",row.names=FALSE,quote=FALSE)

  #Save plots.
  if(PLOTDATA==TRUE){
    ggsave(OUTPUTFILE6,p.filtered.all,width=10,height=20)
  }
}

############################################################################
#SAVE CONTROL FILES:

cat(paste0("############################################################################
Date: ",RUNDATE,"
Job ID: ",JOBID,"
Script: ",SCRIPTNAME,"
Input file: ",INPUTFILE,"
Output file: ",OUTPUTFILE1,"
Output file: ",OUTPUTFILE2,"
Output file: ",OUTPUTFILE3,"
Output file: ",OUTPUTFILE4,"

"),file=paste0(OUTPUTLOCATION,"/README.txt"),append=TRUE)

sink()