############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

plot_cistrans<-function(SNP_CHR=c(1:10), #Vector of SNP chromosomes.
                        SNP_POS=rep(5,10), #Vector of positions of each SNP.
                        SNP_PVALUE=runif(10,min=0,max=0.05), #Vector of p-values for each SNP.
                        TYPE=c(rep("Cis",5),rep("Trans",5)), #Vector of type of eQTL ("Cis","Trans").
                        GENE_ID=sample(paste0("GENE",1:10),10,replace=TRUE), #Vector of gene IDs for each SNP.
                        GENE_CHR=c(1:5,7,8,9,10,6), #Vector of gene chromosomes.
                        GENE_STARTPOS=rep(1,10), #Vector of gene start positions.
                        GENE_ENDPOS=rep(10,10), #Vector of gene end positions.
                        CHRSET){ #Vector of chromosome set.)
  
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
  
  #Format data.
  DATA<-data.frame(SNP_ID=paste0(SNP_CHR,":",SNP_POS),
                   SNP_CHR=SNP_CHR,
                   SNP_POS=SNP_POS,
                   SNP_PVALUE=SNP_PVALUE,
                   TYPE=TYPE,
                   GENE_ID=GENE_ID,
                   GENE_CHR=GENE_CHR,
                   GENE_STARTPOS=GENE_STARTPOS,
                   GENE_ENDPOS=GENE_ENDPOS)
  
  #Define target chromosome set.
  if(missing(CHRSET)){
    CHRSET<-mixedsort(unique(c(DATA$SNP_CHR,DATA$GENE_CHR)))
  }
  DATA<-DATA[DATA$SNP_CHR %in% CHRSET & DATA$GENE_CHR %in% CHRSET,]
  
  #Calculate number of genes per SNP and number of SNPs per gene.
  DATA<-DATA %>%
    group_by(SNP_ID,TYPE) %>%
    mutate(SNP_NGENES=length(unique(GENE_ID)))
  DATA<-DATA %>%
    group_by(GENE_ID,TYPE) %>%
    mutate(GENE_NSNPS=length(unique(SNP_ID))) %>%
    as.data.frame
  #DATA<-DATA %>%
  #  group_by(GENE_ID) %>%
  #  mutate(GENE_NCISSNPS=length(unique(SNP_ID[TYPE=="Cis"]))) %>%
  #  mutate(GENE_NTRANSSNPS=length(unique(SNP_ID[TYPE=="Trans"]))) %>%
  #  ungroup %>%
  #  as.data.frame
  
  #Add missing chromosomes to data set.
  NCHR<-length(CHRSET)
  for(C in CHRSET){
    if(sum(DATA$SNP_CHR==C)<2 || length(unique(DATA$SNP_ID[DATA$SNP_CHR==C]))<2 || 
       sum(DATA$GENE_CHR==C)<2 || length(unique(DATA$GENE_ID[DATA$GENE_CHR==C]))<2){
      TMP<-data.frame(SNP_ID=paste0(C,":",c(50,150)),
                      SNP_CHR=C,
                      SNP_POS=c(50,150),
                      SNP_PVALUE=1,
                      TYPE="DUMMYSNP",
                      SNP_NGENES=0,
                      GENE_ID=paste0("DUMMYGENE-",C,":",c("1:100","101:200")),
                      GENE_CHR=C,
                      GENE_STARTPOS=c(1,101),
                      GENE_ENDPOS=c(100,200),
                      GENE_NSNPS=0)#,
                      #GENE_NCISSNPS=0,
                      #GENE_NTRANSSNPS=0
                      #)
      DATA<-rbind(DATA,TMP)
    }
  }
  
  #Add gene mid position.
  DATA<-DATA %>%
    mutate(GENE_MIDPOS=(DATA$GENE_ENDPOS+DATA$GENE_STARTPOS)/2)
  
  #Add intra-chromosome SNP distance to gene center.
  DATA<-DATA %>%
    mutate(SNPDISTGENE=ifelse(SNP_CHR==GENE_CHR,
                              DATA$SNP_POS-DATA$GENE_MIDPOS,
                              NA))
  
  #Sort by SNP chromosome and position and add plot SNP axis position.
  STDCHRLENGTH<-100
  DATA<-DATA %>%
    arrange(SNP_CHR,SNP_POS)
  SNPxPOS<-data.frame(SNP_ID=unique(DATA$SNP_ID),
                      xPosition=NA)
  rownames(SNPxPOS)<-SNPxPOS$SNP_ID
  for(C in 1:length(CHRSET)){
    CHR_SNPNAMES<-DATA %>%
      filter(SNP_CHR==CHRSET[C]) %>%
      dplyr::select(SNP_ID) %>%
      distinct %>%
      arrange
    SNPxPOS[CHR_SNPNAMES$SNP_ID,2]<-seq(from=(C-1)*STDCHRLENGTH+1,
                                        to=C*STDCHRLENGTH,
                                        length.out=nrow(CHR_SNPNAMES))
  }
  DATA$SNPAXIS_POS<-SNPxPOS[DATA$SNP_ID,2]
  
  #Sort by gene chromosome and position and add plot gene axis position.
  DATA<-DATA %>%
    arrange(GENE_CHR,GENE_MIDPOS)
  GENExPOS<-data.frame(GENE_ID=unique(DATA$GENE_ID),
                       xPosition=NA)
  rownames(GENExPOS)<-GENExPOS$GENE_ID
  for(C in 1:length(CHRSET)){
    CHR_GENENAMES<-DATA %>%
      filter(GENE_CHR==CHRSET[C]) %>%
      dplyr::select(GENE_ID) %>%
      distinct 
    GENExPOS[CHR_GENENAMES$GENE_ID,2]<-seq(from=(C-1)*STDCHRLENGTH+1,
                                           to=C*STDCHRLENGTH,
                                           length.out=nrow(CHR_GENENAMES))
  }
  DATA$GENEAXIS_POS<-GENExPOS[DATA$GENE_ID,2]
  
  #Calculate range of gray/white background.
  RECTCOLORS<-rep(c("gray","white"),NCHR)
  SNP_RECTX<-data.frame(Chromosome=CHRSET,
                        Rectxmin=NA,
                        Rectxmax=NA,
                        Rectcolor=NA)
  GENE_RECTY<-data.frame(Chromosome=CHRSET,
                         Rectymin=NA,
                         Rectymax=NA,
                         Rectcolor=NA)
  for(C in 1:NCHR){
    SNP_RECTX$Rectxmin[C]<-min(DATA$SNPAXIS_POS[DATA$SNP_CHR==SNP_RECTX$Chromosome[C]])
    SNP_RECTX$Rectxmax[C]<-max(DATA$SNPAXIS_POS[DATA$SNP_CHR==SNP_RECTX$Chromosome[C]])
    SNP_RECTX$Rectcolor[C]<-RECTCOLORS[C]
    GENE_RECTY$Rectymin[C]<-min(DATA$GENEAXIS_POS[DATA$GENE_CHR==GENE_RECTY$Chromosome[C]])
    GENE_RECTY$Rectymax[C]<-max(DATA$GENEAXIS_POS[DATA$GENE_CHR==GENE_RECTY$Chromosome[C]])
    GENE_RECTY$Rectcolor[C]<-RECTCOLORS[C]
  }

  #Calculate center of each chromosome on the x axis.
  SNP_AXIS_SET<-DATA %>% 
    group_by(SNP_CHR) %>%
    summarize(center=(max(SNPAXIS_POS)+min(SNPAXIS_POS))/2)
  GENE_AXIS_SET<-DATA %>% 
    group_by(GENE_CHR) %>%
    summarize(center=(max(GENEAXIS_POS)+min(GENEAXIS_POS))/2)
  
  ##########################################
  #Create plots.
  #Define chromosome colors.
  CHRCOLORS<-rep(c("deepskyblue","darkblue"),NCHR)
  CISCHRCOLORS<-CHRCOLORS[match(DATA$SNP_CHR,CHRSET)]
  CHRCOLORS<-rep(c("tomato","darkred"),NCHR)
  TRANSCHRCOLORS<-CHRCOLORS[match(DATA$SNP_CHR,CHRSET)]
  
  if(("Cis" %in% unique(DATA$TYPE)) & ("Trans" %in% unique(DATA$TYPE))){
    
    #Manhattan plots.
    #Cis-eQTLs.
    p1a<-ggplot(data=DATA[DATA$TYPE=="Cis",])+
      geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
                fill=SNP_RECTX$Rectcolor,
                alpha=0.3)+
      geom_point(aes(x=SNPAXIS_POS,y=-log2(SNP_PVALUE)),
                 color=CISCHRCOLORS[DATA$TYPE=="Cis"],
                 shape=18,
                 size=0.5)+
      scale_x_continuous(label=SNP_AXIS_SET$SNP_CHR,breaks=SNP_AXIS_SET$center)+
      scale_y_continuous(limits=c(min(-log2(DATA$SNP_PVALUE[DATA$TYPE=="Cis"]),na.rm=TRUE),
                                  max(-log2(DATA$SNP_PVALUE[DATA$TYPE=="Cis"]),na.rm=TRUE)))+
      labs(x=expression("Chromosome"),y=expression('-Log'[2]~'P-value'))+
      #ggtitle(expression("Cis-SNPs"))+
      annotate("text",
               x=0.9*max(SNP_RECTX$Rectxmax),
               y=0.9*max(-log2(DATA$SNP_PVALUE[DATA$TYPE=="Cis"])),
               label="Cis-SNPs")+
      #facet_wrap(~TYPE)+
      myggplottheme_blank
    #p1a
    
    #Trans-eQTLs.
    p1b<-ggplot(data=DATA[DATA$TYPE=="Trans",])+
      geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
                fill=SNP_RECTX$Rectcolor,
                alpha=0.3)+
      geom_point(aes(x=SNPAXIS_POS,y=-log2(SNP_PVALUE)),
                 color=TRANSCHRCOLORS[DATA$TYPE=="Trans"],
                 shape=18,
                 size=0.5)+
      scale_x_continuous(label=SNP_AXIS_SET$SNP_CHR,breaks=SNP_AXIS_SET$center)+
      scale_y_continuous(limits=c(min(-log2(DATA$SNP_PVALUE[DATA$TYPE=="Trans"]),na.rm=TRUE),
                                  max(-log2(DATA$SNP_PVALUE[DATA$TYPE=="Trans"]),na.rm=TRUE)))+
      labs(x=expression("Chromosome"),y=expression('-Log'[2]~'P-value'))+
      #ggtitle(expression("Trans-SNPs"))+
      annotate("text",
               x=0.9*max(SNP_RECTX$Rectxmax),
               y=0.9*max(-log2(DATA$SNP_PVALUE[DATA$TYPE=="Trans"])),
               label="Trans-SNPs")+
      #facet_wrap(~TYPE)+
      myggplottheme_blank
    #p1b
    
    p1<-ggarrange(p1a,p1b,
                  ncol=1)
    #p1
    
    #Histograms of p-values.
    #Cis-eQTLs.
    HISTCISDATA<-data.frame(Pvalue=DATA$SNP_PVALUE[DATA$TYPE=="Cis"],
                            TYPE="Cis")
    p2a<-ggplot(HISTCISDATA,aes(x=Pvalue))+
      geom_histogram(fill="blue",color="#e9ecef",alpha=0.9)+
      #geom_histogram(aes(y=..density..),fill="blue",color="#e9ecef",alpha=0.9)+
      #geom_density()+
      #geom_vline(aes(xintercept=mean(Pvalue)),
      #           color="black",linetype="dashed",linewidth=0.5)+
      labs(x=expression("P-value (FDR)"),y=expression("# Cis-SNPs"))+
      #ggtitle("Cis-eQTLs")+
      #facet_wrap(~TYPE)+
      myggplottheme+
      theme(axis.text.x=element_text(angle=0,vjust=0,hjust=0.5))
    #p2a
    
    #Trans-eQTLs.
    HISTTRANSDATA<-data.frame(Pvalue=DATA$SNP_PVALUE[DATA$TYPE=="Trans"],
                              TYPE="Trans")
    p2b<-ggplot(HISTTRANSDATA,aes(x=Pvalue))+
      geom_histogram(fill="red",color="#e9ecef",alpha=0.9)+
      #geom_histogram(aes(y=..density..),fill="red",color="#e9ecef",alpha=0.9)+
      #geom_density()+
      #geom_vline(aes(xintercept=mean(Pvalue)),
      #           color="black",linetype="dashed",linewidth=0.5)+
      labs(x=expression("P-value (FDR)"),y=expression("# Trans-SNPs"))+
      #ggtitle("Trans-eQTLs")+
      #facet_wrap(~TYPE)+
      myggplottheme+
      theme(axis.text.x=element_text(angle=0,vjust=0,hjust=0.5))
    #p2b
    
    p2<-ggarrange(p2a,p2b,
                  ncol=1)
    #p2
    
    #Cis-trans plot.
    p3<-ggplot(data=DATA)+
      geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
                fill=SNP_RECTX$Rectcolor,
                alpha=0.3)+
      geom_rect(data=GENE_RECTY,aes(xmin=-Inf,xmax=Inf,ymin=Rectymin,ymax=Rectymax),
                fill=GENE_RECTY$Rectcolor,
                alpha=0.3)+
      geom_point(data=DATA[DATA$TYPE=="Trans",],aes(x=SNPAXIS_POS,y=GENEAXIS_POS),
                 color="red",
                 shape=18,
                 size=0.5)+
      geom_point(data=DATA[DATA$TYPE=="Cis",],aes(x=SNPAXIS_POS,y=GENEAXIS_POS),
                 color="blue",
                 shape=18,
                 size=0.5)+
      geom_point(data=DATA[DATA$TYPE=="DUMMYSNP",],aes(x=SNPAXIS_POS,y=GENEAXIS_POS),
                 color=NA,
                 shape=18,
                 size=2)+
      scale_x_continuous(label=SNP_AXIS_SET$SNP_CHR,breaks=SNP_AXIS_SET$center)+
      scale_y_continuous(label=GENE_AXIS_SET$GENE_CHR,breaks=GENE_AXIS_SET$center)+
      labs(x=expression("SNP position"),y=expression("Gene position"))+
      ggtitle(expression("Cis-trans plot"))+
      myggplottheme_blank
    #p3
    
    #Number of cis-genes per SNP.
    p4a<-unique(DATA[,c("SNP_ID","TYPE","SNPAXIS_POS","SNP_NGENES")]) %>%
      filter(TYPE=="Cis") %>%
      ggplot()+
      geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
                fill=SNP_RECTX$Rectcolor,
                alpha=0.3)+
      geom_bar(aes(x=SNPAXIS_POS,y=SNP_NGENES),
               stat="identity",
               col="blue")+
      geom_hline(aes(yintercept=mean(SNP_NGENES)),
                 color="black",linetype="dashed",linewidth=0.5)+
      scale_x_continuous(label=SNP_AXIS_SET$SNP_CHR,breaks=SNP_AXIS_SET$center)+
      labs(x=expression("SNP position"),y=expression("# Cis-genes"))+
      #ggtitle("Number of genes per SNP")+
      myggplottheme_blank+
      theme(axis.title.x=element_blank(),
            axis.text.y=element_text(angle=90,vjust=1,hjust=0.5))
    #p4a
    
    #Number of trans-genes per SNP.
    p4b<-unique(DATA[,c("SNP_ID","TYPE","SNPAXIS_POS","SNP_NGENES")]) %>%
      filter(TYPE=="Trans") %>%
      ggplot()+
      geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
                fill=SNP_RECTX$Rectcolor,
                alpha=0.3)+
      geom_bar(aes(x=SNPAXIS_POS,y=SNP_NGENES),
               stat="identity",
               col="red")+
      geom_hline(aes(yintercept=mean(SNP_NGENES)),
                 color="black",linetype="dashed",linewidth=0.5)+
      scale_x_continuous(label=SNP_AXIS_SET$SNP_CHR,breaks=SNP_AXIS_SET$center)+
      labs(x=expression("SNP position"),y=expression("# Trans-genes"))+
      #ggtitle("Number of genes per SNP")+
      myggplottheme_blank+
      theme(axis.title.x=element_blank(),
            axis.text.y=element_text(angle=90,vjust=1,hjust=0.5))
    #p4b
    
    p4<-ggarrange(p4a,p4b,
                  ncol=1)
    #p4
    
    ##Number of SNPs per gene.
    #p8<-unique(DATA[,c("GENE_ID","GENEAXIS_POS","GENE_NSNPS")]) %>%
    #  ggplot()+
    #  geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
    #            fill=SNP_RECTX$Rectcolor,
    #            alpha=0.3)+
    #  geom_bar(aes(x=GENEAXIS_POS,y=GENE_NSNPS),
    #           stat="identity",
    #           col="orange")+
    #  geom_hline(aes(yintercept=mean(SNP_NSNPS)),
    #             color="black",linetype="dashed",linewidth=0.5)+
    #  scale_x_continuous(label=GENE_AXIS_SET$GENE_CHR,breaks=GENE_AXIS_SET$center)+
    #  labs(x="Chromosome",y="# SNPs")+
    #  #ggtitle("Number of SNPs per gene")+
    #  myggplottheme_blank+
    #  theme(axis.title.y=element_blank(),
    #        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5))+
    #  rotate()
    ##p8
    
    #Number of cis-SNPs per gene.
    p5a<-unique(DATA[,c("GENE_ID","TYPE","GENEAXIS_POS","GENE_NSNPS")]) %>%
      filter(TYPE=="Cis") %>%
      ggplot()+
      geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
                fill=SNP_RECTX$Rectcolor,
                alpha=0.3)+
      geom_bar(aes(x=GENEAXIS_POS,y=GENE_NSNPS),
               stat="identity",
               col="blue")+
      geom_hline(aes(yintercept=mean(GENE_NSNPS)),
                 color="black",linetype="dashed",linewidth=0.5)+
      scale_x_continuous(label=GENE_AXIS_SET$GENE_CHR,breaks=GENE_AXIS_SET$center)+
      labs(x=expression("Chromosome"),y=expression("# Cis-SNPs"))+
      #ggtitle("Number of cis-SNPs per gene")+
      myggplottheme_blank+
      theme(axis.title.y=element_blank(),
            axis.text.x=element_text(angle=0,vjust=1,hjust=0.5)
      )+
      rotate()
    #p5a
    
    #Number of trans-SNPs per gene.
    p5b<-unique(DATA[,c("GENE_ID","TYPE","GENEAXIS_POS","GENE_NSNPS")]) %>%
      filter(TYPE=="Trans") %>%
      ggplot()+
      geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
                fill=SNP_RECTX$Rectcolor,
                alpha=0.3)+
      geom_bar(aes(x=GENEAXIS_POS,y=GENE_NSNPS),
               stat="identity",
               col="red")+
      geom_hline(aes(yintercept=mean(GENE_NSNPS)),
                 color="black",linetype="dashed",linewidth=0.5)+
      scale_x_continuous(label=GENE_AXIS_SET$GENE_CHR,breaks=GENE_AXIS_SET$center)+
      labs(x=expression("Chromosome"),y=expression("# Trans-SNPs"))+
      #ggtitle("Number of trans-SNPs per gene")+
      myggplottheme_blank+
      theme(axis.title.y=element_blank(),
            axis.text.x=element_text(angle=0,vjust=1,hjust=0.5)
      )+
      rotate()
    #p5b
    
    p5<-ggarrange(p5a,p5b,
                  nrow=1)
    #p5
    
    #Intra-chromosomal distance between SNP and gene center.
    p6<-ggplot(data=DATA)+
    geom_point(data=DATA[DATA$SNP_CHR %in% CHRSET &
                             DATA$GENE_CHR %in% CHRSET &
                             DATA$TYPE=="Trans",],aes(x=SNPDISTGENE,y=-log2(SNP_PVALUE)),
                 color="red",
                 shape=18,
                 size=1)+
      geom_point(data=DATA[DATA$SNP_CHR %in% CHRSET &
                             DATA$GENE_CHR %in% CHRSET &
                             DATA$TYPE=="Cis",],aes(x=SNPDISTGENE,y=-log2(SNP_PVALUE)),
                 color="blue",
                 shape=18,
                 size=0.5)+
      labs(x=expression("Distance betwen SNP and gene center (bp)"),y=expression('-Log'[2]~'P-value'))+
      #ggtitle(expression("Distance between SNP and gene center"))+
      myggplottheme_blank+
      theme(axis.text.x=element_text(angle=0,vjust=0,hjust=0.5))
    #p6
    
    #Correlation between number of cis- and trans-eQTLs per gene.
    NSNPSxNGENESDATA<-unique(DATA[,c("SNP_ID","GENE_ID","TYPE")]) %>%
      group_by(GENE_ID) %>%
      mutate(GENE_NCISSNPS=length(unique(SNP_ID[TYPE=="Cis"]))) %>%
      mutate(GENE_NTRANSSNPS=length(unique(SNP_ID[TYPE=="Trans"]))) %>%
      ungroup %>%
      as.data.frame
    p7<-ggplot(data=unique(NSNPSxNGENESDATA[,c("GENE_ID","GENE_NCISSNPS","GENE_NTRANSSNPS")]))+
      geom_point(aes(x=GENE_NCISSNPS,y=GENE_NTRANSSNPS,group=GENE_ID),
                 color="gray30",
                 shape=18,
                 size=1)+
      labs(x=expression("# Cis-SNPs"),y=expression("# Trans-SNPs"))+
      #ggtitle("Cis- vs. trans-SNPs per gene")+
      myggplottheme_blank+
      theme(axis.text.x=element_text(angle=0,vjust=0,hjust=0.5))
    #p7
    
    p.all<-ggarrange(ggarrange(p1,p2,p6,
                               nrow=3,ncol=1,
                               heights=c(2,2,2),
                               labels=c("A","B","C")),
                     ggarrange(ggarrange(p3,p4,
                                         nrow=2,
                                         heights=c(4,2),
                                         labels=c("D","E")),
                               ggarrange(p5,p7,
                                         nrow=2,
                                         heights=c(4,2),
                                         labels=c("F","G")),
                               ncol=2,nrow=1,align="v",
                               heights=c(4,2),
                               widths=c(4,2)),
                     ncol=2,nrow=1,
                     widths=c(1,2.5))
    #p.all
    
  }else{
    
    #Manhattan plots.
    p1<-ggplot(data=DATA)+
      geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
                fill=SNP_RECTX$Rectcolor,
                alpha=0.3)+
      geom_point(aes(x=SNPAXIS_POS,y=-log2(SNP_PVALUE)),
                 color=CISCHRCOLORS,
                 shape=18,
                 size=0.5)+
      scale_x_continuous(label=SNP_AXIS_SET$SNP_CHR,breaks=SNP_AXIS_SET$center)+
      scale_y_continuous(limits=c(min(-log2(DATA$SNP_PVALUE),na.rm=TRUE),
                                  max(-log2(DATA$SNP_PVALUE),na.rm=TRUE)))+
      labs(x=expression("Chromosome"),y=expression('-Log'[2]~'P-value'))+
      #ggtitle(expression("Cis-SNPs"))+
      annotate("text",
               x=0.9*max(SNP_RECTX$Rectxmax),
               y=0.9*max(-log2(DATA$SNP_PVALUE)),
               label="All SNPs")+
      #facet_wrap(~TYPE)+
      myggplottheme_blank
    #p1
    
    #Histograms of p-values.
    HISTCISDATA<-data.frame(Pvalue=DATA$SNP_PVALUE)
    p2<-ggplot(HISTCISDATA,aes(x=Pvalue))+
      geom_histogram(fill="blue",color="#e9ecef",alpha=0.9)+
      #geom_histogram(aes(y=..density..),fill="blue",color="#e9ecef",alpha=0.9)+
      #geom_density()+
      #geom_vline(aes(xintercept=mean(Pvalue)),
      #           color="black",linetype="dashed",linewidth=0.5)+
      labs(x=expression("P-value (FDR)"),y=expression("# SNPs"))+
      #ggtitle("Cis-eQTLs")+
      #facet_wrap(~TYPE)+
      myggplottheme+
      theme(axis.text.x=element_text(angle=0,vjust=0,hjust=0.5))
    #p2
    
    #Cis-trans plot.
    p3<-ggplot(data=DATA)+
      geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
                fill=SNP_RECTX$Rectcolor,
                alpha=0.3)+
      geom_rect(data=GENE_RECTY,aes(xmin=-Inf,xmax=Inf,ymin=Rectymin,ymax=Rectymax),
                fill=GENE_RECTY$Rectcolor,
                alpha=0.3)+
      geom_point(data=DATA[DATA$TYPE!="DUMMYSNP",],aes(x=SNPAXIS_POS,y=GENEAXIS_POS),
                 color="blue",
                 shape=18,
                 size=0.5)+
      geom_point(data=DATA[DATA$TYPE=="DUMMYSNP",],aes(x=SNPAXIS_POS,y=GENEAXIS_POS),
                 color=NA,
                 shape=18,
                 size=2)+
      scale_x_continuous(label=SNP_AXIS_SET$SNP_CHR,breaks=SNP_AXIS_SET$center)+
      scale_y_continuous(label=GENE_AXIS_SET$GENE_CHR,breaks=GENE_AXIS_SET$center)+
      labs(x=expression("SNP position"),y=expression("Gene position"))+
      ggtitle(expression("Cis-trans plot"))+
      myggplottheme_blank
    #p3
    
    #Number of genes per SNP.
    p4<-unique(DATA[,c("SNP_ID","SNPAXIS_POS","SNP_NGENES")]) %>%
      ggplot()+
      geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
                fill=SNP_RECTX$Rectcolor,
                alpha=0.3)+
      geom_bar(aes(x=SNPAXIS_POS,y=SNP_NGENES),
               stat="identity",
               col="blue")+
      geom_hline(aes(yintercept=mean(SNP_NGENES)),
                 color="black",linetype="dashed",linewidth=0.5)+
      scale_x_continuous(label=SNP_AXIS_SET$SNP_CHR,breaks=SNP_AXIS_SET$center)+
      labs(x=expression("SNP position"),y=expression("# Genes"))+
      #ggtitle("Number of genes per SNP")+
      myggplottheme_blank+
      theme(axis.title.x=element_blank(),
            axis.text.y=element_text(angle=90,vjust=1,hjust=0.5))
    #p4
    
    ##Number of SNPs per gene.
    #p8<-unique(DATA[,c("GENE_ID","GENEAXIS_POS","GENE_NSNPS")]) %>%
    #  ggplot()+
    #  geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
    #            fill=SNP_RECTX$Rectcolor,
    #            alpha=0.3)+
    #  geom_bar(aes(x=GENEAXIS_POS,y=GENE_NSNPS),
    #           stat="identity",
    #           col="orange")+
    #  geom_hline(aes(yintercept=mean(SNP_NSNPS)),
    #             color="black",linetype="dashed",linewidth=0.5)+
    #  scale_x_continuous(label=GENE_AXIS_SET$GENE_CHR,breaks=GENE_AXIS_SET$center)+
    #  labs(x="Chromosome",y="# SNPs")+
    #  #ggtitle("Number of SNPs per gene")+
    #  myggplottheme_blank+
    #  theme(axis.title.y=element_blank(),
    #        axis.text.x=element_text(angle=0,vjust=1,hjust=0.5))+
    #  rotate()
    ##p8
    
    #Number of SNPs per gene.
    p5<-unique(DATA[,c("GENE_ID","GENEAXIS_POS","GENE_NSNPS")]) %>%
      ggplot()+
      geom_rect(data=SNP_RECTX,aes(xmin=Rectxmin,xmax=Rectxmax,ymin=-Inf,ymax=Inf),
                fill=SNP_RECTX$Rectcolor,
                alpha=0.3)+
      geom_bar(aes(x=GENEAXIS_POS,y=GENE_NSNPS),
               stat="identity",
               col="blue")+
      geom_hline(aes(yintercept=mean(GENE_NSNPS)),
                 color="black",linetype="dashed",linewidth=0.5)+
      scale_x_continuous(label=GENE_AXIS_SET$GENE_CHR,breaks=GENE_AXIS_SET$center)+
      labs(x=expression("Chromosome"),y=expression("# SNPs"))+
      #ggtitle("Number of SNPs per gene")+
      myggplottheme_blank+
      theme(axis.title.y=element_blank(),
            axis.text.x=element_text(angle=0,vjust=1,hjust=0.5)
      )+
      rotate()
    #p5
    
    #Intra-chromosomal distance between SNP and gene center.
    p6<-ggplot(data=DATA)+
      geom_point(data=DATA[DATA$SNP_CHR %in% CHRSET &
                             DATA$GENE_CHR %in% CHRSET,],aes(x=SNPDISTGENE,y=-log2(SNP_PVALUE)),
                 color="blue",
                 shape=18,
                 size=0.5)+
      labs(x=expression("Distance betwen SNP and gene center (bp)"),y=expression('-Log'[2]~'P-value'))+
      #ggtitle(expression("Distance between SNP and gene center"))+
      myggplottheme_blank+
      theme(axis.text.x=element_text(angle=0,vjust=0,hjust=0.5))
    #p6
    
    p.all<-ggarrange(ggarrange(p1,p2,p6,
                               nrow=3,ncol=1,
                               heights=c(2,2,2),
                               labels=c("A","B","C")),
                     ggarrange(ggarrange(p3,p4,
                                         nrow=2,
                                         heights=c(4,2),
                                         labels=c("D","E")),
                               ggarrange(p5,
                                         nrow=2,
                                         heights=c(4,2),
                                         labels=c("F",NULL)),
                               ncol=2,nrow=1,align="v",
                               heights=c(4,2),
                               widths=c(4,2)),
                     ncol=2,nrow=1,
                     widths=c(1,2.5))
    #p.all

  }
  
  return(p.all)
  
}


  
























