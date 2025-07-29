############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

plot_genotype_phenotype<-function(MATRIX_GENOTYPES, #Data frame of genotypes with SNPs per row (1st column is SNP ID) and samples per column.
                                  MATRIX_PHENOTYPES, #Data frame of phenotypic values with phenotype per row (1st column is phenotype ID) and samples per column.
                                  VECTOR_SNP_ID, #Vector of SNPs to be plotted (up to 25).
                                  VECTOR_PHENOTYPE_ID, #Vector of phenotypes to be plotted (up to 25).
                                  MINGENFREQ=0, #Minimum genotype frequency (<1) or count (>1).
                                  PLOTTITLE="Genotype x Phenotype",
                                  XLABEL="Genotype",
                                  YLABEL="Phenotype",
                                  MAXNPLOTS=25){ #Maximum number of plots.
  
  #Load libraries.
  library(ggplot2)
  library(tidyverse)
  library(gtools)
  
  myggplottheme_blank_nolegend<-theme(title=element_text(size=10,face="bold"),
                                      axis.title=element_text(size=10,face="bold"),
                                      axis.text=element_text(size=10),
                                      axis.text.x=element_text(angle=60,size=8,vjust=0.5),
                                      legend.position="none",
                                      #legend.title=element_text(size=10,face="bold"),
                                      #legend.text=element_text(size=10),
                                      #legend.key=element_blank(),
                                      panel.grid=element_line(colour="gray90"),
                                      panel.grid.major.x=element_blank(),
                                      panel.grid.minor.x=element_blank(),
                                      panel.background=element_rect(fill="white",colour="black"),
                                      panel.grid.major=element_blank(),
                                      panel.grid.minor=element_blank(),
                                      strip.background=element_rect(colour="black",
                                                                    fill="white"))
  
  #Subset genotype data.
  DATA_GENOTYPES<-MATRIX_GENOTYPES[MATRIX_GENOTYPES[,1] %in% VECTOR_SNP_ID,]
  
  #Subset phenotype data.
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
  
  #Data frame with target SNP-phenotype pairs.
  SNP_PHENOTYPE<-data.frame(SNP_ID=VECTOR_SNP_ID,
                            PHENOTYPE_ID=VECTOR_PHENOTYPE_ID)
  SNP_PHENOTYPE<-SNP_PHENOTYPE %>%
    filter(!(SNP_ID %in% SNPltMINGENFREQ))
  
  if(nrow(SNP_PHENOTYPE)>MAXNPLOTS){
    SNP_PHENOTYPE<-SNP_PHENOTYPE[1:min(MAXNPLOTS,nrow(SNP_PHENOTYPE)),]
  }
  
  #Create data frame with genotypes and phenotypes.
  DATA_GENOTYPE_PHENOTYPE<-data.frame(SNP=NA,
                                      PHENOTYPE=NA,
                                      SAMPLE=NA,
                                      GENOTYPE=NA,
                                      PHENOTYPEVALUE=NA)
  
  for(i in 1:nrow(SNP_PHENOTYPE)){
    TMP<-data.frame(SNP=SNP_PHENOTYPE$SNP_ID[i],
                    PHENOTYPE=SNP_PHENOTYPE$PHENOTYPE_ID[i],
                    SAMPLE=colnames(DATA_GENOTYPES)[2:ncol(DATA_GENOTYPES)],
                    GENOTYPE=unname(unlist(DATA_GENOTYPES[
                      DATA_GENOTYPES[,1]==SNP_PHENOTYPE$SNP_ID[i],
                      2:ncol(DATA_GENOTYPES)])),
                    PHENOTYPEVALUE=as.numeric(unname(unlist(DATA_PHENOTYPES[
                      DATA_PHENOTYPES[,1]==SNP_PHENOTYPE$PHENOTYPE_ID[i],
                      2:ncol(DATA_PHENOTYPES)]))))
    DATA_GENOTYPE_PHENOTYPE<-rbind(DATA_GENOTYPE_PHENOTYPE,TMP)
  }
  DATA_GENOTYPE_PHENOTYPE<-DATA_GENOTYPE_PHENOTYPE[!is.na(DATA_GENOTYPE_PHENOTYPE$SNP),]
  
  #Create plot.
  p.genxphen<-DATA_GENOTYPE_PHENOTYPE %>%
    ggplot(aes(x=GENOTYPE,y=PHENOTYPEVALUE))+ 
    geom_jitter(color="black",shape=20,size=0.2,
                height=0,width=0.3)+
    geom_boxplot(aes(fill=GENOTYPE),
                 width=0.5,color="black",alpha=0.7,
                 outlier.color="black",outlier.shape=20,outlier.size=0.2)+
    labs(x=XLABEL,y=YLABEL)+
    facet_wrap(~SNP+PHENOTYPE,
               nrow=floor(sqrt(nrow(SNP_PHENOTYPE))),
               scales="free")+
    ggtitle(PLOTTITLE)+
    myggplottheme_blank_nolegend
  p.genxphen
  
}
