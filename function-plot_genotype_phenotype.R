############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

plot_genotype_phenotype<-function(MATRIX_GENOTYPES, #Data frame of genotypes with SNPs per row (1st column is SNP ID) and samples per column.
                                  MATRIX_PHENOTYPES, #Data frame of trait values with trait per row (1st column is trait ID) and samples per column.
                                  VECTOR_SNP_ID, #Vector of SNPs to be plotted (up to 25).
                                  VECTOR_PHENOTYPE_ID, #Vector of traits to be plotted (up to 25).
                                  MINGENFREQ=0, #Minimum genotype frequency (<1) or count (>1).
                                  RESCALEOUTLIERS=FALSE, #Replace outlier values with Q1-1.5*IQR or Q3+1.5*IQR.
                                  MAXNPLOTS=25, #Maximum number of plots.
                                  PLOTTITLE="Genotype x Phenotype",
                                  XLABEL="Genotype",
                                  YLABEL="Phenotype",
                                  COLOR_BOX="gray40",
                                  COLOR_OUTLIER="red"){ 
  
  #Load libraries.
  library(ggplot2)
  library(tidyverse)
  library(gtools)
  
  myggplottheme_blank_nolegend<-theme(title=element_text(size=10,face="bold"),
                                      axis.title=element_text(size=10,face="bold"),
                                      axis.text=element_text(size=10),
                                      axis.text.x=element_text(angle=60,size=8,vjust=0.5),
                                      legend.position="none",
                                      panel.grid=element_line(colour="gray90"),
                                      panel.grid.major.x=element_blank(),
                                      panel.grid.minor.x=element_blank(),
                                      panel.background=element_rect(fill="white",colour="black"),
                                      panel.grid.major=element_blank(),
                                      panel.grid.minor=element_blank(),
                                      strip.background=element_rect(colour="black",
                                                                    fill="white"))
  
  #Subset genotype data. Keep only SNPs that are in target SNP vector.
  DATA_GENOTYPES<-MATRIX_GENOTYPES[MATRIX_GENOTYPES[,1] %in% VECTOR_SNP_ID,]
  
  #Subset phenotype data. Keep only traits that are in target trait vector.
  DATA_PHENOTYPES<-MATRIX_PHENOTYPES[MATRIX_PHENOTYPES[,1] %in% VECTOR_PHENOTYPE_ID,]
  
  #Transform phased genotypes into unphased genotypes.
  DATA_GENOTYPES<-data.frame(lapply(DATA_GENOTYPES,function(x){
    x<-gsub("\\|","/",x)
    x<-gsub("C/A","A/C",x)
    x<-gsub("G/A","A/G",x)
    x<-gsub("T/A","A/T",x)
    x<-gsub("G/C","C/G",x)
    x<-gsub("T/C","C/T",x)
    gsub("T/G","G/T",x)}))
  rownames(DATA_GENOTYPES)<-DATA_GENOTYPES[,1]
  
  #Get common samples.
  GENOTYPES_SAMPLES<-colnames(DATA_GENOTYPES)[2:ncol(DATA_GENOTYPES)]
  PHENOTYPES_SAMPLES<-colnames(DATA_PHENOTYPES)[2:ncol(DATA_PHENOTYPES)]
  SAMPLES_COMMON<-GENOTYPES_SAMPLES[GENOTYPES_SAMPLES %in% PHENOTYPES_SAMPLES]
  
  #Match order of columns in genotype and phenotype data.
  DATA_GENOTYPES<-DATA_GENOTYPES[,c(colnames(DATA_GENOTYPES)[1],SAMPLES_COMMON)]
  DATA_PHENOTYPES<-DATA_PHENOTYPES[,c(colnames(DATA_PHENOTYPES)[1],SAMPLES_COMMON)]
  
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
  if(MINGENFREQ<1){
    GENOTYPES_COUNTS$N_Genotypes_ltMINGENFREQ<-rowSums(GENOTYPES_COUNTS[,2:12]>0 & (GENOTYPES_COUNTS[,2:12]/GENOTYPES_COUNTS$N_Samples)<MINGENFREQ)
    #SNPs that have at least one genotype with frequency <MINGENFREQ.
    SNPltMINGENFREQ<-GENOTYPES_COUNTS$SNP[which(rowSums(GENOTYPES_COUNTS[,2:12]>0 & (GENOTYPES_COUNTS[,2:12]/GENOTYPES_COUNTS$N_Samples)<MINGENFREQ)>0)]
  }else{
    GENOTYPES_COUNTS$N_Genotypes_ltMINGENFREQ<-rowSums(GENOTYPES_COUNTS[,2:12]>0 & GENOTYPES_COUNTS[,2:12]<MINGENFREQ)
    #SNPs that have at least one genotype with frequency <MINGENFREQ.
    SNPltMINGENFREQ<-GENOTYPES_COUNTS$SNP[which(rowSums(GENOTYPES_COUNTS[,2:12]>0 & GENOTYPES_COUNTS[,2:12]<MINGENFREQ)>0)]
  }
  
  #Filter out SNPs that have at least one genotype with frequency <MINGENFREQ.
  DATA_GENOTYPES<-DATA_GENOTYPES %>%
    filter(!(DATA_GENOTYPES[,1] %in% SNPltMINGENFREQ))
  if(nrow(DATA_GENOTYPES)==0){return(paste0("Genotypes do not meet the MINGENFREQ=",MINGENFREQ," requirement."))}
  
  #Data frame with target SNP-phenotype pairs.
  SNP_PHENOTYPE<-data.frame(SNP_ID=VECTOR_SNP_ID,
                            PHENOTYPE_ID=VECTOR_PHENOTYPE_ID)
  SNP_PHENOTYPE<-SNP_PHENOTYPE %>%
    filter(!(SNP_ID %in% SNPltMINGENFREQ))
  
  if(nrow(SNP_PHENOTYPE)>MAXNPLOTS){
    SNP_PHENOTYPE<-SNP_PHENOTYPE[1:min(MAXNPLOTS,nrow(SNP_PHENOTYPE)),]
  }
  
  #Create data frame with genotypes and phenotypes.
  DATA_GENOTYPE_PHENOTYPE<-data.frame(SNP_ID=NA,
                                      Phenotype=NA,
                                      Sample=NA,
                                      Genotype=NA,
                                      PhenotypeValue=NA)
  
  for(i in 1:nrow(SNP_PHENOTYPE)){
    TMP<-data.frame(SNP_ID=SNP_PHENOTYPE$SNP_ID[i],
                    Phenotype=SNP_PHENOTYPE$PHENOTYPE_ID[i],
                    Sample=colnames(DATA_GENOTYPES)[2:ncol(DATA_GENOTYPES)],
                    Genotype=unname(unlist(DATA_GENOTYPES[
                      DATA_GENOTYPES[,1]==SNP_PHENOTYPE$SNP_ID[i],
                      2:ncol(DATA_GENOTYPES)])),
                    PhenotypeValue=as.numeric(unname(unlist(DATA_PHENOTYPES[
                      DATA_PHENOTYPES[,1]==SNP_PHENOTYPE$PHENOTYPE_ID[i],
                      2:ncol(DATA_PHENOTYPES)]))))
    DATA_GENOTYPE_PHENOTYPE<-rbind(DATA_GENOTYPE_PHENOTYPE,TMP)
  }
  #Remove NAs.
  DATA_GENOTYPE_PHENOTYPE<-DATA_GENOTYPE_PHENOTYPE[!is.na(DATA_GENOTYPE_PHENOTYPE$SNP_ID),]
  
  #Find outliers per SNP, trait and genotype.
  DATA_GENOTYPE_PHENOTYPE<-DATA_GENOTYPE_PHENOTYPE %>%
    group_by(SNP_ID,Phenotype,Genotype) %>%
    mutate(Q1=quantile(PhenotypeValue,0.25,na.rm=TRUE),
           Q3=quantile(PhenotypeValue,0.75,na.rm=TRUE),
           IQR=Q3-Q1,
           Lower=Q1-1.5*IQR,
           Upper=Q3+1.5*IQR,
           Outlier=(PhenotypeValue<Lower | PhenotypeValue>Upper),
           PhenotypeValue_RescaledOutlier=case_when(PhenotypeValue<Lower~Lower,
                                                    PhenotypeValue>Upper~Upper,
                                                    TRUE~PhenotypeValue)) %>%
    ungroup()
  
  #Create plot.
  if(RESCALEOUTLIERS==TRUE){
    p.genxphen<-DATA_GENOTYPE_PHENOTYPE %>%
      ggplot(aes(x=Genotype,y=PhenotypeValue_RescaledOutlier))+ 
      geom_boxplot(aes(fill=Genotype),
                   width=0.5,color=COLOR_BOX,alpha=0.7,
                   outlier.color="black",outlier.shape=NA,outlier.size=0.2)+
      #geom_jitter(color="black",shape=20,size=0.2,
      #            height=0,width=0.3)+
      geom_jitter(aes(y=PhenotypeValue_RescaledOutlier,color=Outlier,shape=Outlier),
                  size=0.4,height=0,width=0.3)+
      scale_color_manual(values=c("FALSE"="black","TRUE"=COLOR_OUTLIER))+
      scale_shape_manual(values=c("FALSE"=20,"TRUE"=17))+
      labs(x=XLABEL,y=YLABEL)+
      facet_wrap(~SNP_ID+Phenotype,
                 nrow=floor(sqrt(nrow(SNP_PHENOTYPE))),
                 scales="free",
                 drop=TRUE)+
      ggtitle(PLOTTITLE)+
      myggplottheme_blank_nolegend
  }else{
    p.genxphen<-DATA_GENOTYPE_PHENOTYPE %>%
      ggplot(aes(x=Genotype,y=PhenotypeValue))+ 
      geom_boxplot(aes(fill=Genotype),
                   width=0.5,color=COLOR_BOX,alpha=0.7,
                   outlier.color="black",outlier.shape=NA,outlier.size=0.2)+
      #geom_jitter(color="black",shape=20,size=0.2,
      #            height=0,width=0.3)+
      geom_jitter(aes(y=PhenotypeValue,color=Outlier,shape=Outlier),
                  size=0.4,height=0,width=0.3)+
      scale_color_manual(values=c("FALSE"="black","TRUE"=COLOR_OUTLIER))+
      scale_shape_manual(values=c("FALSE"=20,"TRUE"=17))+
      labs(x=XLABEL,y=YLABEL)+
      facet_wrap(~SNP_ID+Phenotype,
                 nrow=floor(sqrt(nrow(SNP_PHENOTYPE))),
                 scales="free",
                 drop=TRUE)+
      ggtitle(PLOTTITLE)+
      myggplottheme_blank_nolegend
  }
  
  p.genxphen
  
}

plot_genotype_phenotype_lmresiduals<-function(MATRIX_GENOTYPES, #Data frame of genotypes with SNPs per row (1st column is SNP ID) and samples per column.
                                              MATRIX_PHENOTYPES, #Data frame of trait values with trait per row (1st column is trait ID) and samples per column.
                                              VECTOR_SNP_ID, #Vector of SNPs to be plotted (up to 25). VECTOR_SNP_ID and VECTOR_PHENOTYPE_ID must have the same length.
                                              VECTOR_PHENOTYPE_ID, #Vector of traits to be plotted (up to 25).
                                              VECTOR_COVARIATES, #Vector of coma-separated covariates for each trait in VECTOR_PHENOTYPE_ID. Covariates must be comma-separated and all covariates must be present in MATRIX_PHENOTYPES. 
                                              MINGENFREQ=0, #Minimum genotype frequency (<1) or count (>1).
                                              RESCALEOUTLIERS=FALSE, #Replace outlier values with Q1-1.5*IQR or Q3+1.5*IQR.
                                              RESIDUALS_STANDARDIZED=TRUE,
                                              MAXNPLOTS=25, #Maximum number of plots.
                                              OUTPUTRESIDUALS=FALSE,
                                              PLOTTITLE="Genotype x Standardized LM residuals",
                                              XLABEL="Genotype",
                                              YLABEL="Standardized LM residuals",
                                              COLOR_BOX="gray40",
                                              COLOR_OUTLIER="red"){ 
  
  #Load libraries.
  library(ggplot2)
  library(tidyverse)
  library(gtools)
  
  myggplottheme_blank_nolegend<-theme(title=element_text(size=10,face="bold"),
                                      axis.title=element_text(size=10,face="bold"),
                                      axis.text=element_text(size=10),
                                      axis.text.x=element_text(angle=60,size=8,vjust=0.5),
                                      legend.position="none",
                                      panel.grid=element_line(colour="gray90"),
                                      panel.grid.major.x=element_blank(),
                                      panel.grid.minor.x=element_blank(),
                                      panel.background=element_rect(fill="white",colour="black"),
                                      panel.grid.major=element_blank(),
                                      panel.grid.minor=element_blank(),
                                      strip.background=element_rect(colour="black",
                                                                    fill="white"))
  
  #Subset genotype data. Keep only SNPs that are in target SNP vector.
  DATA_GENOTYPES<-MATRIX_GENOTYPES[MATRIX_GENOTYPES[,1] %in% VECTOR_SNP_ID,]
  
  #Subset phenotype data. Keep only traits that are in target trait vector.
  LIST_PHENOTYPES_COVARIATES<-unique(c(VECTOR_PHENOTYPE_ID,
                                       unique(trimws(unlist(strsplit(VECTOR_COVARIATES,","))))))
  DATA_PHENOTYPES<-MATRIX_PHENOTYPES[MATRIX_PHENOTYPES[,1] %in% LIST_PHENOTYPES_COVARIATES,]
  
  #Transform phased genotypes into unphased genotypes.
  DATA_GENOTYPES<-data.frame(lapply(DATA_GENOTYPES,function(x){
    x<-gsub("\\|","/",x)
    x<-gsub("C/A","A/C",x)
    x<-gsub("G/A","A/G",x)
    x<-gsub("T/A","A/T",x)
    x<-gsub("G/C","C/G",x)
    x<-gsub("T/C","C/T",x)
    gsub("T/G","G/T",x)}))
  rownames(DATA_GENOTYPES)<-DATA_GENOTYPES[,1]
  
  #Get common samples.
  GENOTYPES_SAMPLES<-colnames(DATA_GENOTYPES)[2:ncol(DATA_GENOTYPES)]
  PHENOTYPES_SAMPLES<-colnames(DATA_PHENOTYPES)[2:ncol(DATA_PHENOTYPES)]
  SAMPLES_COMMON<-GENOTYPES_SAMPLES[GENOTYPES_SAMPLES %in% PHENOTYPES_SAMPLES]
  
  #Match order of columns in genotype and phenotype data.
  DATA_GENOTYPES<-DATA_GENOTYPES[,c(colnames(DATA_GENOTYPES)[1],SAMPLES_COMMON)]
  DATA_PHENOTYPES<-DATA_PHENOTYPES[,c(colnames(DATA_PHENOTYPES)[1],SAMPLES_COMMON)]
  
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
  if(MINGENFREQ<1){
    GENOTYPES_COUNTS$N_Genotypes_ltMINGENFREQ<-rowSums(GENOTYPES_COUNTS[,2:12]>0 & (GENOTYPES_COUNTS[,2:12]/GENOTYPES_COUNTS$N_Samples)<MINGENFREQ)
    #SNPs that have at least one genotype with frequency <MINGENFREQ.
    SNPltMINGENFREQ<-GENOTYPES_COUNTS$SNP[which(rowSums(GENOTYPES_COUNTS[,2:12]>0 & (GENOTYPES_COUNTS[,2:12]/GENOTYPES_COUNTS$N_Samples)<MINGENFREQ)>0)]
  }else{
    GENOTYPES_COUNTS$N_Genotypes_ltMINGENFREQ<-rowSums(GENOTYPES_COUNTS[,2:12]>0 & GENOTYPES_COUNTS[,2:12]<MINGENFREQ)
    #SNPs that have at least one genotype with frequency <MINGENFREQ.
    SNPltMINGENFREQ<-GENOTYPES_COUNTS$SNP[which(rowSums(GENOTYPES_COUNTS[,2:12]>0 & GENOTYPES_COUNTS[,2:12]<MINGENFREQ)>0)]
  }
  
  #Filter out SNPs that have at least one genotype with frequency <MINGENFREQ.
  DATA_GENOTYPES<-DATA_GENOTYPES %>%
    filter(!(DATA_GENOTYPES[,1] %in% SNPltMINGENFREQ))
  if(nrow(DATA_GENOTYPES)==0){return(paste0("Genotypes do not meet the MINGENFREQ=",MINGENFREQ," requirement."))}
  
  #Data frame with target SNP-phenotype pairs.
  SNP_PHENOTYPE<-data.frame(SNP_ID=VECTOR_SNP_ID,
                            PHENOTYPE_ID=VECTOR_PHENOTYPE_ID,
                            COVARIATES=VECTOR_COVARIATES)
  SNP_PHENOTYPE<-SNP_PHENOTYPE %>%
    filter(!(SNP_ID %in% SNPltMINGENFREQ))
  
  #Set maximum number of plots.
  if(nrow(SNP_PHENOTYPE)>MAXNPLOTS){
    SNP_PHENOTYPE<-SNP_PHENOTYPE[1:min(MAXNPLOTS,nrow(SNP_PHENOTYPE)),]
  }
  
  DATA_PHENOTYPES_RESIDUALS<-DATA_PHENOTYPES[0,]
  for(i in 1:nrow(SNP_PHENOTYPE)){
    
    SNP<-SNP_PHENOTYPE$SNP_ID[i]
    PHENOTYPE<-SNP_PHENOTYPE$PHENOTYPE_ID[i]
    
    #Get phenotype-specific list of covariates.
    LMCOV<-trimws(strsplit(SNP_PHENOTYPE$COVARIATES[i],",")[[1]])
    
    #Get data.
    LMDATA<-data.frame(Sample=SAMPLES_COMMON,
                       Phenotype=unname(unlist(DATA_PHENOTYPES[DATA_PHENOTYPES[,1]==PHENOTYPE,SAMPLES_COMMON])),
                       Genotype=unname(unlist(DATA_GENOTYPES[DATA_GENOTYPES[,1]==SNP,SAMPLES_COMMON])))
    LMDATA<-cbind(LMDATA,t(DATA_PHENOTYPES[LMCOV,SAMPLELIST]))
    rownames(LMDATA)<-LMDATA$Sample
    
    #Define numeric variables.
    for(c in c(2,4:ncol(LMDATA))){
      LMDATA[,c]<-as.numeric(LMDATA[,c])
    }
    
    #Linear model.
    OUT_LM<-lm(LMDATA[,2]~.,data=LMDATA[,-c(1,2)])
    
    #Add residuals to data frame.
    DATA_PHENOTYPES_RESIDUALS[i,]<-NA
    DATA_PHENOTYPES_RESIDUALS$SNP_ID[i]<-SNP
    DATA_PHENOTYPES_RESIDUALS$PHENOTYPE_ID[i]<-PHENOTYPE
    if(RESIDUALS_STANDARDIZED==TRUE){
      DATA_PHENOTYPES_RESIDUALS[i,names(rstandard(OUT_LM))]<-rstandard(OUT_LM)
    }else{
      DATA_PHENOTYPES_RESIDUALS[i,names(OUT_LM$residuals)]<-OUT_LM$residuals
    }
    DATA_PHENOTYPES_RESIDUALS$LinearModel[i]<-paste0(PHENOTYPE," ~ ",paste0(c("Genotype",LMCOV),collapse=" + "))
    
  }
  
  #Reorder columns.
  TMP<-cbind(DATA_PHENOTYPES_RESIDUALS[,c("PHENOTYPE_ID","SNP_ID","LinearModel")],DATA_PHENOTYPES_RESIDUALS[,2:(ncol(DATA_PHENOTYPES_RESIDUALS)-2)])
  DATA_PHENOTYPES_RESIDUALS<-TMP
  rm(TMP)
  
  #Make residuals numeric.
  for(c in 4:ncol(DATA_PHENOTYPES_RESIDUALS)){DATA_PHENOTYPES_RESIDUALS[,c]<-as.numeric(DATA_PHENOTYPES_RESIDUALS[,c])}
  
  #Create data frame with genotypes and phenotype residuals.
  DATA_GENOTYPE_PHENOTYPE<-DATA_PHENOTYPES_RESIDUALS %>%
    pivot_longer(cols=4:ncol(DATA_PHENOTYPES_RESIDUALS),
                 names_to="Sample",
                 values_to="PhenotypeValue")
  
  #Add genotype column.
  DATA_GENOTYPE_PHENOTYPE$Genotype<-NA
  for(SNP in unique(DATA_GENOTYPE_PHENOTYPE$SNP_ID)){
    LIST_SAMPLES<-DATA_GENOTYPE_PHENOTYPE$Sample[DATA_GENOTYPE_PHENOTYPE$SNP_ID==SNP]
    LIST_GENOTYPES<-unname(unlist(DATA_GENOTYPES[DATA_GENOTYPES[,1]==SNP,LIST_SAMPLES]))
    DATA_GENOTYPE_PHENOTYPE$Genotype[DATA_GENOTYPE_PHENOTYPE$SNP_ID==SNP]<-LIST_GENOTYPES
  }
  
  #Remove NAs.
  DATA_GENOTYPE_PHENOTYPE<-DATA_GENOTYPE_PHENOTYPE[!is.na(DATA_GENOTYPE_PHENOTYPE$Genotype),]
  DATA_GENOTYPE_PHENOTYPE<-DATA_GENOTYPE_PHENOTYPE[!is.na(DATA_GENOTYPE_PHENOTYPE$PhenotypeValue),]
  
  #Find outliers per SNP, trait and genotype.
  DATA_GENOTYPE_PHENOTYPE<-DATA_GENOTYPE_PHENOTYPE %>%
    group_by(SNP_ID,PHENOTYPE_ID,Genotype) %>%
    mutate(Q1=quantile(PhenotypeValue,0.25,na.rm=TRUE),
           Q3=quantile(PhenotypeValue,0.75,na.rm=TRUE),
           IQR=Q3-Q1,
           Lower=Q1-1.5*IQR,
           Upper=Q3+1.5*IQR,
           Outlier=(PhenotypeValue<Lower | PhenotypeValue>Upper),
           PhenotypeValue_RescaledOutlier=case_when(PhenotypeValue<Lower~Lower,
                                                    PhenotypeValue>Upper~Upper,
                                                    TRUE~PhenotypeValue)) %>%
    ungroup()
  
  #Create plot.
  if(RESCALEOUTLIERS==TRUE){
    p.genxphen<-DATA_GENOTYPE_PHENOTYPE %>%
      ggplot(aes(x=Genotype,y=PhenotypeValue_RescaledOutlier))+ 
      geom_boxplot(aes(fill=Genotype),
                   width=0.5,color=COLOR_BOX,alpha=0.7,
                   outlier.color="black",outlier.shape=NA,outlier.size=0.2)+
      #geom_jitter(color="black",shape=20,size=0.2,
      #            height=0,width=0.3)+
      geom_jitter(aes(y=PhenotypeValue_RescaledOutlier,color=Outlier,shape=Outlier),
                  size=0.4,height=0,width=0.3)+
      scale_color_manual(values=c("FALSE"="black","TRUE"=COLOR_OUTLIER))+
      scale_shape_manual(values=c("FALSE"=20,"TRUE"=17))+
      labs(x=XLABEL,y=YLABEL)+
      facet_wrap(~SNP_ID+PHENOTYPE_ID,
                 nrow=floor(sqrt(nrow(SNP_PHENOTYPE))),
                 scales="free",
                 drop=TRUE)+
      ggtitle(PLOTTITLE)+
      myggplottheme_blank_nolegend
  }else{
    p.genxphen<-DATA_GENOTYPE_PHENOTYPE %>%
      ggplot(aes(x=Genotype,y=PhenotypeValue))+ 
      geom_boxplot(aes(fill=Genotype),
                   width=0.5,color=COLOR_BOX,alpha=0.7,
                   outlier.color="black",outlier.shape=NA,outlier.size=0.2)+
      #geom_jitter(color="black",shape=20,size=0.2,
      #            height=0,width=0.3)+
      geom_jitter(aes(y=PhenotypeValue,color=Outlier,shape=Outlier),
                  size=0.4,height=0,width=0.3)+
      scale_color_manual(values=c("FALSE"="black","TRUE"=COLOR_OUTLIER))+
      scale_shape_manual(values=c("FALSE"=20,"TRUE"=17))+
      labs(x=XLABEL,y=YLABEL)+
      facet_wrap(~SNP_ID+PHENOTYPE_ID,
                 nrow=floor(sqrt(nrow(SNP_PHENOTYPE))),
                 scales="free",
                 drop=TRUE)+
      ggtitle(PLOTTITLE)+
      myggplottheme_blank_nolegend
  }
  
  if(OUTPUTRESIDUALS==TRUE){
    return(list(DATA_RESIDUALS=DATA_GENOTYPE_PHENOTYPE,
                PLOT=p.genxphen))
  }else{
    return(p.genxphen)
  }
  
}

############################################################################
############################### DEPRECATED #################################
############################################################################
plot_genotype_phenotype_pair<-function(MATRIX_GENOTYPES, #Data frame of genotypes with SNPs per row (1st column is SNP ID) and samples per column.
                                       MATRIX_PHENOTYPES, #Data frame of phenotypic values with phenotype per row (1st column is phenotype ID, 2nd column is SNP ID) and samples per column.
                                       MINGENFREQ=0, #Minimum genotype frequency (<1) or count (>1).
                                       PLOTTITLE="Genotype x Phenotype",
                                       XLABEL="Genotype",
                                       YLABEL="Phenotype",
                                       MAXNPLOTS=25, #Maximum number of plots.
                                       BOXCOLOR="black"){ 
  #This function should be used when phenotypic values are SNP-specific (e.g., phenotypic values are residuals from a linear model involving a specific SNP).
  
  #Load libraries.
  library(ggplot2)
  library(tidyverse)
  library(gtools)
  
  myggplottheme_blank_nolegend<-theme(title=element_text(size=10,face="bold"),
                                      axis.title=element_text(size=10,face="bold"),
                                      axis.text=element_text(size=10),
                                      axis.text.x=element_text(angle=60,size=8,vjust=0.5),
                                      legend.position="none",
                                      panel.grid=element_line(colour="gray90"),
                                      panel.grid.major.x=element_blank(),
                                      panel.grid.minor.x=element_blank(),
                                      panel.background=element_rect(fill="white",colour="black"),
                                      panel.grid.major=element_blank(),
                                      panel.grid.minor=element_blank(),
                                      strip.background=element_rect(colour="black",
                                                                    fill="white"))
  
  colnames(MATRIX_PHENOTYPES)[1:2]<-c("PHENOTYPE_ID","SNP_ID")
  colnames(MATRIX_GENOTYPES)[1]<-"SNP_ID"
  
  #Subset genotype data.
  DATA_GENOTYPES<-MATRIX_GENOTYPES[MATRIX_GENOTYPES[,1] %in% MATRIX_PHENOTYPES[,2],]
  
  #Subset phenotype data.
  DATA_PHENOTYPES<-MATRIX_PHENOTYPES[MATRIX_PHENOTYPES[,2] %in% DATA_GENOTYPES[,1],]
  
  #Transform phased genotypes into unphased genotypes.
  DATA_GENOTYPES<-data.frame(lapply(DATA_GENOTYPES,function(x){
    x<-gsub("\\|","/",x)
    x<-gsub("C/A","A/C",x)
    x<-gsub("G/A","A/G",x)
    x<-gsub("T/A","A/T",x)
    x<-gsub("G/C","C/G",x)
    x<-gsub("T/C","C/T",x)
    gsub("T/G","G/T",x)}))
  rownames(DATA_GENOTYPES)<-DATA_GENOTYPES[,1]
  
  #Get common samples.
  GENOTYPES_SAMPLES<-colnames(DATA_GENOTYPES)[2:ncol(DATA_GENOTYPES)]
  PHENOTYPES_SAMPLES<-colnames(DATA_PHENOTYPES)[3:ncol(DATA_PHENOTYPES)]
  SAMPLES_COMMON<-GENOTYPES_SAMPLES[GENOTYPES_SAMPLES %in% PHENOTYPES_SAMPLES]
  
  #Match order of columns in genotype and phenotype data.
  DATA_GENOTYPES<-DATA_GENOTYPES[,c(colnames(DATA_GENOTYPES)[1],SAMPLES_COMMON)]
  DATA_PHENOTYPES<-DATA_PHENOTYPES[,c(colnames(DATA_PHENOTYPES)[1:2],SAMPLES_COMMON)]
  
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
                               TT=rowSums(DATA_GENOTYPES=="T/T"))
  GENOTYPES_COUNTS$N_Samples<-rowSums(GENOTYPES_COUNTS[,2:11])
  GENOTYPES_COUNTS$N_Genotypes<-rowSums(GENOTYPES_COUNTS[,2:11]>0)
  if(MINGENFREQ<1){
    GENOTYPES_COUNTS$N_Genotypes_ltMINGENFREQ<-rowSums(GENOTYPES_COUNTS[,2:11]>0 & (GENOTYPES_COUNTS[,2:11]/GENOTYPES_COUNTS$N_Samples)<MINGENFREQ)
    #SNPs that have at least one genotype with frequency <MINGENFREQ.
    SNPltMINGENFREQ<-GENOTYPES_COUNTS$SNP[which(rowSums(GENOTYPES_COUNTS[,2:11]>0 & (GENOTYPES_COUNTS[,2:11]/GENOTYPES_COUNTS$N_Samples)<MINGENFREQ)>0)]
  }else{
    GENOTYPES_COUNTS$N_Genotypes_ltMINGENFREQ<-rowSums(GENOTYPES_COUNTS[,2:11]>0 & GENOTYPES_COUNTS[,2:11]<MINGENFREQ)
    #SNPs that have at least one genotype with frequency <MINGENFREQ.
    SNPltMINGENFREQ<-GENOTYPES_COUNTS$SNP[which(rowSums(GENOTYPES_COUNTS[,2:11]>0 & GENOTYPES_COUNTS[,2:11]<MINGENFREQ)>0)]
  }

  #Filter out SNPs that have at least one genotype with frequency <MINGENFREQ.
  DATA_GENOTYPES<-DATA_GENOTYPES %>%
    filter(!(DATA_GENOTYPES[,1] %in% SNPltMINGENFREQ))
  if(nrow(DATA_GENOTYPES)==0){return(paste0("Genotypes do not meet the MINGENFREQ=",MINGENFREQ," requirement."))}
  
  DATA_PHENOTYPES<-DATA_PHENOTYPES %>%
    filter(!(DATA_PHENOTYPES[,2] %in% SNPltMINGENFREQ))
  if(nrow(DATA_PHENOTYPES)==0){return(paste0("Genotypes do not meet the MINGENFREQ=",MINGENFREQ," requirement."))}
  
  #Set maximum number of plots.
  if(nrow(DATA_PHENOTYPES)>MAXNPLOTS){
    DATA_PHENOTYPES<-DATA_PHENOTYPES[1:min(MAXNPLOTS,nrow(DATA_PHENOTYPES)),]
  }
  
  #Create data frame with genotypes and phenotypes.
  DATA_GENOTYPE_PHENOTYPE<-data.frame(SNP_ID=NA,
                                      Phenotype=NA,
                                      Sample=NA,
                                      Genotype=NA,
                                      PhenotypeValue=NA)
  
  for(i in 1:nrow(DATA_PHENOTYPES)){
    TMP<-data.frame(SNP_ID=DATA_PHENOTYPES$SNP_ID[i],
                    Phenotype=DATA_PHENOTYPES$PHENOTYPE_ID[i],
                    Sample=colnames(DATA_GENOTYPES)[2:ncol(DATA_GENOTYPES)],
                    Genotype=unname(unlist(DATA_GENOTYPES[
                      DATA_GENOTYPES[,1]==DATA_PHENOTYPES$SNP_ID[i],
                      2:ncol(DATA_GENOTYPES)])),
                    PhenotypeValue=as.numeric(unname(unlist(DATA_PHENOTYPES[i,3:ncol(DATA_PHENOTYPES)]))))
    DATA_GENOTYPE_PHENOTYPE<-rbind(DATA_GENOTYPE_PHENOTYPE,TMP)
  }
  DATA_GENOTYPE_PHENOTYPE<-DATA_GENOTYPE_PHENOTYPE[!is.na(DATA_GENOTYPE_PHENOTYPE$SNP_ID),]
  
  #Create plot.
  p.genxphen<-DATA_GENOTYPE_PHENOTYPE %>%
    ggplot(aes(x=Genotype,y=PhenotypeValue))+ 
    geom_boxplot(aes(fill=Genotype),
                 width=0.5,color=BOXCOLOR,alpha=0.7,
                 outlier.color="black",outlier.shape=NA,outlier.size=0.2)+
    geom_jitter(color="black",shape=20,size=0.2,
                height=0,width=0.3)+
    labs(x=XLABEL,y=YLABEL)+
    facet_wrap(~SNP_ID+Phenotype,
               nrow=floor(sqrt(nrow(DATA_PHENOTYPES))),
               scales="free")+
    ggtitle(PLOTTITLE)+
    myggplottheme_blank_nolegend
  p.genxphen
  
}
