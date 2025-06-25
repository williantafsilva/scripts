############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

#Manhattan plot.
plot_manhattan<-function(
  CHRVECTOR=paste0("Chr",1:10), #Vector of chromosomes.
  POSVECTOR=1:10, #Vector of positions.
  CHRSET, #Set of chromosomes to show in the plots.
  VALUES=runif(10), #Vector of y values (p-values, likelihood values, iHS values, etc).
  REGIONS, #Data frame (Chromosome,RegionStart,RegionEnd) with chromosomal regions to be shaded in per chromosome plots.
  ABOVETHRESHOLD=Inf, #Threshold above which points are plotted with alpha=1. If both ABOVETHRESHOLD and BELOWTHRESHOLD are set to Inf and -Inf, all points are plotted with alpha=1.
  BELOWTHRESHOLD=-Inf, #Threshold below which points are plotted with alpha=1.
  COLORS=c("deepskyblue","darkblue"), #Alternating colors for adjacent chromosomes.
  STDCHRLENGTH=TRUE, #Use standardized chromosome length?
  PLOTPERCHR=FALSE, #Plot per chromosome? (This options never uses standardized chromosome length).
  XLABEL="Chromosome", #X-axis label.
  YLABEL="Value", #Y-axis label.
  PLOTTITLE="Manhattan plot"){ #
  
  #Load libraries.
  library(ggplot2)
  library(tidyverse)
  library(gtools)
  
  #Plot theme.
  myggplottheme<-theme(title=element_text(size=10,face="bold"),
                       axis.title=element_text(size=10,face="bold"),
                       axis.text=element_text(size=10),
                       axis.text.x=element_text(angle=60,size=8,vjust=0.5),
                       legend.position="none",
                       legend.title=element_text(size=10,face="bold"),
                       legend.text=element_text(size=10),
                       legend.key=element_blank(),
                       panel.grid=element_line(colour="gray90"),
                       panel.grid.major.x=element_blank(),
                       panel.grid.minor.x=element_blank(),
                       panel.background=element_rect(fill="white",colour="black"),
                       panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                       strip.background=element_rect(colour="black",
                                                     fill="white"))
  
  #Format data.
  DATA<-data.frame(CHR=CHRVECTOR,
                   POS=POSVECTOR,
                   VALUE=VALUES)

  #Define target chromosome set.
  if(missing(CHRSET)){
    CHRSET<-mixedsort(unique(as.character(DATA$CHR)))
  }
  DATA<-DATA[DATA$CHR %in% CHRSET,]
  
  #Calculate average value per chromosome.
  DATA<-DATA %>%
    group_by(CHR) %>%
    mutate(AVERAGE=mean(VALUE))
  
  #Sort data by chromosome name and position.
  DATA<-DATA[with(DATA,order(CHR,POS)),]
  
  #Add missing chromosomes to data set.
  NCHR<-length(CHRSET)
  for(C in CHRSET){
    if(sum(DATA$CHR==C)<2){
      TMP<-data.frame(CHR=C,
                      POS=c(25,75),
                      VALUE=NA,
                      AVERAGE=NA)
      DATA<-rbind(DATA,TMP)
    }
  }
  
  #Sort data by chromosome name and position.
  DATA<-DATA[with(DATA,order(CHR,POS)),]
  
  #Standardize chromosome length (for plotting aesthetics purposes; standardize chromosome length).
  DATA$xPOSITION<-NA
  if(STDCHRLENGTH){
    STDLENGTH<-100
    for(C in 1:NCHR){
      DATA$xPOSITION[DATA$CHR==CHRSET[C]]<-
        seq(from=(C-1)*STDLENGTH+1,
            to=C*STDLENGTH,
            length.out=sum(DATA$CHR==CHRSET[C]))
    }
  }else{
    DATA$xPOSITION<-1:nrow(DATA)
  }
  
  #Define transparency (alpha value).
  if(ABOVETHRESHOLD<Inf | BELOWTHRESHOLD>-Inf){
    DATA$ALPHA<-0.3
  }else{
    DATA$ALPHA<-1
  }
  if(ABOVETHRESHOLD<Inf){
    DATA$ALPHA[DATA$VALUE>=ABOVETHRESHOLD]<-1
  }
  if(BELOWTHRESHOLD>-Inf){
    DATA$ALPHA[DATA$VALUE<=BELOWTHRESHOLD]<-1
  }
  
  #Calculate range of gray/white background.
  RECTCOLORS<-rep(c("gray","white"),NCHR)
  RECTX<-data.frame(CHR=CHRSET,
                    Rectxmin=NA,
                    Rectxmax=NA,
                    Rectcolor=NA)
  for(C in 1:NCHR){
    RECTX$Rectxmin[C]<-min(DATA$xPOSITION[DATA$CHR==RECTX$CHR[C]])
    RECTX$Rectxmax[C]<-max(DATA$xPOSITION[DATA$CHR==RECTX$CHR[C]])
    RECTX$Rectcolor[C]<-RECTCOLORS[C]
  }
  
  #Define chromosome colors.
  CHRCOLORS<-rep(COLORS,NCHR)
  DATA$CHR_COLOR<-CHRCOLORS[match(DATA$CHR,CHRSET)]

  #Calculate center of each chromosome on the x axis.
  AXIS_SET<-DATA %>% 
    group_by(CHR) %>%
    summarize(CHR_CENTER=mean(xPOSITION))
  
  #Plot.
  if(PLOTPERCHR==FALSE){
    p<-ggplot(data=DATA)+
      geom_rect(data=RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
                fill=RECTX$Rectcolor,
                alpha=0.3)+
      geom_point(aes(xPOSITION,VALUE),
                 color=DATA$CHR_COLOR,
                 alpha=DATA$ALPHA,
                 shape=18,
                 size=1)+
      scale_x_continuous(label=AXIS_SET$CHR,breaks=AXIS_SET$CHR_CENTER)+
      scale_y_continuous(limits=c(min(DATA$VALUE,na.rm=TRUE),max(DATA$VALUE,na.rm=TRUE)))+
      geom_hline(yintercept=BELOWTHRESHOLD,color="black",linetype="dashed")+
      geom_hline(yintercept=ABOVETHRESHOLD,color="black",linetype="dashed")+
      labs(x=XLABEL,y=YLABEL)+
      ggtitle(PLOTTITLE)+
      myggplottheme
  }else{
    if(missing(REGIONS)){
      p<-ggplot(data=DATA)+
        geom_point(aes(POS,VALUE),
                   color=DATA$CHR_COLOR,
                   alpha=DATA$ALPHA,
                   shape=18,
                   size=1)+
        #scale_x_continuous(label=AXIS_SET$Chromosome,breaks=AXIS_SET$center)+
        scale_y_continuous(limits=c(min(DATA$VALUE,na.rm=TRUE),max(DATA$VALUE,na.rm=TRUE)))+
        geom_hline(yintercept=BELOWTHRESHOLD,color="black",linetype="dashed")+
        geom_hline(yintercept=ABOVETHRESHOLD,color="black",linetype="dashed")+
        labs(x=XLABEL,y=YLABEL)+
        facet_wrap(~CHR,nrow=ceiling(sqrt(NCHR)),scales="free_x")+
        ggtitle(PLOTTITLE)+
        myggplottheme
    }else{
      colnames(REGIONS)[1:3]<-c("CHR","RegionStart","RegionEnd")
      p<-ggplot(data=DATA)+
        geom_rect(data=REGIONS,aes(xmin=RegionStart,xmax=RegionEnd,ymin=-Inf,ymax=Inf),
                  fill="gray",
                  alpha=0.3,
                  inherit.aes=FALSE)+
        geom_point(aes(POS,VALUE),
                   color=DATA$CHR_COLOR,
                   alpha=DATA$ALPHA,
                   shape=18,
                   size=1)+
        #scale_x_continuous(label=AXIS_SET$Chromosome,breaks=AXIS_SET$center)+
        scale_y_continuous(limits=c(min(DATA$VALUE,na.rm=TRUE),max(DATA$VALUE,na.rm=TRUE)))+
        geom_hline(yintercept=BELOWTHRESHOLD,color="black",linetype="dashed")+
        geom_hline(yintercept=ABOVETHRESHOLD,color="black",linetype="dashed")+
        labs(x=XLABEL,y=YLABEL)+
        facet_wrap(~CHR,nrow=ceiling(sqrt(NCHR)),scales="free_x")+
        ggtitle(PLOTTITLE)+
        myggplottheme
    }
    
  }

  return(p)
  
}
