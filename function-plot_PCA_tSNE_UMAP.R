############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

#Principal components analysis.
plot_PCA_tSNE_UMAP<-function(inputdata=data.frame(Feature=c("A","B","C","D","E"),
                                                  Sample1=runif(5),
                                                  Sample2=runif(5),
                                                  Sample3=runif(5),
                                                  Sample4=runif(5),
                                                  row.names=1),
                             inputmetadata=data.frame(Group=c("Sex","Age"),
                                                      Sample1=c("M","Juvenile"),
                                                      Sample2=c("F","Adult"),
                                                      Sample3=c("M","Adult"),
                                                      Sample4=c("M","Juvenile"),
                                                      row.names=1),
                             colorgroup="Age",
                             method="PCA"){
  
  #Load libraries.
  library(tidyverse)
  library(ggfortify)
  library(ggplot2)
  library(PCAtools)
  library(Rtsne)
  library(umap)
  library(factoextra)
  library(corrplot)
  library(cluster)
  library(dplyr)
  library(readr)
  library(tidyr)

  #Plot theme.
  myggplottheme<-theme(title=element_text(size=10,face="bold"),
                       axis.title=element_text(size=10,face="bold"),
                       axis.text=element_text(size=10),
                       axis.text.x=element_text(angle=60,size=8,vjust=0.5),
                       #legend.position="none",
                       legend.title=element_text(size=10,face="bold"),
                       legend.text=element_text(size=10),
                       legend.key=element_blank(),
                       panel.grid=element_line(colour="gray90"),
                       #panel.grid.major.x=element_blank(),
                       #panel.grid.minor.x=element_blank(),
                       panel.background=element_rect(fill="white",colour="black"),
                       #panel.grid.major=element_blank(),
                       #panel.grid.minor=element_blank(),
                       strip.background=element_rect(colour="black",
                                                     fill="white"))
  
  #Load data.
  #DATA<-read.table(inputdata,header=TRUE,sep="\t") %>% column_to_rownames(var=colnames(.)[1])
  DATA<-inputdata
  DATA<-t(DATA) #Transpose.
  DATA<-DATA[,which(apply(DATA,2,var)!=0)] #Remove zero-variance features.

  #METADATA<-read.table(inputmetadata,header=TRUE,sep="\t") %>% column_to_rownames(var=colnames(.)[1])
  METADATA<-inputmetadata
  METADATA[is.na(METADATA)]<-"Unkown"
  METADATA<-t(METADATA)
  
  NCOLORS<-length(unique(DATA_PLOT[,colorgroup]))
  
  if(method=="PCA"){
    #Perform PCA.
    DATA_PCA<-prcomp(DATA[,which(apply(DATA,2,var)!=0)],scale.=TRUE)
    DATA_PC1_PC2<-as.data.frame(DATA_PCA$x[,1:2])
    colnames(DATA_PC1_PC2)<-c("PC1","PC2")
    #Create plots.
    DATA_PLOT<-cbind(DATA_PC1_PC2,METADATA)
    PLOT<-ggplot(DATA_PLOT,aes(x=PC1,y=PC2,color=.data[[colorgroup]]))+
      geom_point()+
      scale_colour_manual(values=rainbow(NCOLORS))+
      myggplottheme
  }
  
  if(method=="tSNE"){
    #Perform t-SNE.
    PERPLEXITY<-30
    DATA_tSNE<-Rtsne(DATA,perplexity=PERPLEXITY,check_duplicates=FALSE)
    DATA_tSNE1_tSNE2<-as.data.frame(DATA_tSNE$Y)
    colnames(DATA_tSNE1_tSNE2)<-c("tSNE1","tSNE2")
    #Create plots.
    DATA_PLOT<-cbind(DATA_tSNE1_tSNE2,METADATA)
    PLOT<-ggplot(DATA_PLOT,aes(x=tSNE1,y=tSNE2,color=.data[[colorgroup]]))+
      geom_point()+
      scale_colour_manual(values=rainbow(NCOLORS))+
      myggplottheme
  }

  if(method=="UMAP"){
    #Perform UMAP.
    DATA_UMAP<-umap(DATA)
    DATA_UMAP1_UMAP2<-as.data.frame(DATA_UMAP$layout)
    colnames(DATA_UMAP1_UMAP2)<-c("UMAP1","UMAP2")
    #Create plots.
    DATA_PLOT<-cbind(DATA_UMAP1_UMAP2,METADATA)
    PLOT<-ggplot(DATA_PLOT,aes(x=UMAP1,y=UMAP2,color=.data[[colorgroup]]))+
      geom_point()+
      scale_colour_manual(values=rainbow(NCOLORS))+
      myggplottheme
  }
  
  return(PLOT)
  
}
