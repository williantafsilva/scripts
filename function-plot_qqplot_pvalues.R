############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

#Manhattan plot.
plot_qqplot_pvalues<-function(
    PVALUES=c(runif(1000,min=0,max=1),runif(500,min=0,max=0.01)), #Vector of p values.
    CONFINTERVAL=0.95, #Confidence level of the point-wise confidence band. The default is 0.95. The confidence intervals are computed under the assumption of the p-values being drawn independently from a uniform [0,1] distribution. 
    SIGTHRESHOLD=0.01, #Significance level.
    LOG10SCALE=TRUE, #Plot log10 scale.
    PLOTTITLE="QQ plot of p-values",
    MAXNPOINTS=Inf){ #Maximum number of points to plot (sampled randomly).
  
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
  
  N<-length(PVALUES)
  
  #Calculate lambda statistic (measure of inflation). 
  #The lambda statistic should be close to 1 if the points fall within the expected range, 
  #or greater than one if the observed p-values are more significant than expected.
  #If analysis results your data follows the normal chi-squared distribution 
  #(no inflation), the expected λ value is 1. If the λ value is greater than 1, then 
  #this may be evidence for some systematic bias that needs to be corrected in your analysis.
  CHISQ<-qchisq(1-PVALUES,1)
  LAMBDA<-median(CHISQ)/qchisq(0.5,1)
  
  if(LOG10SCALE==TRUE){
    
    #Format data.
    DATA<-data.frame(Observed=-log10(sort(PVALUES)),
                     Expected=-log10(ppoints(N)),
                     CI_lower=-log10(qbeta(p=(1-CONFINTERVAL)/2,shape1=1:N,shape2=N:1)),
                     CI_upper=-log10(qbeta(p=(1+CONFINTERVAL)/2,shape1=1:N,shape2=N:1)))
    if(MAXNPOINTS<Inf & nrow(DATA)>MAXNPOINTS){ #Sample points to be plotted.
      DATA<-DATA[sample(1:nrow(DATA),MAXNPOINTS,replace=FALSE),]
    }
    
    p<-ggplot(DATA)+
      geom_ribbon(mapping=aes(x=Expected,ymin=CI_lower,ymax=CI_upper),
                  alpha=0.1)+
      geom_segment(aes(x=0,y=-log10(SIGTHRESHOLD),xend=-log10(SIGTHRESHOLD),yend=-log10(SIGTHRESHOLD)),
                   color="red",linetype=3,linewidth=0.5,alpha=0.5)+
      geom_segment(aes(x=-log10(SIGTHRESHOLD),y=0,xend=-log10(SIGTHRESHOLD),yend=-log10(SIGTHRESHOLD)),
                   color="red",linetype=3,linewidth=0.5,alpha=0.5)+
      geom_point(aes(Expected,Observed),shape=16,size=1)+
      geom_abline(intercept=0,slope=1,color="blue",alpha=0.5)+
      geom_line(aes(Expected,CI_upper),linetype=3,linewidth=0.5)+
      geom_line(aes(Expected,CI_lower),linetype=3,linewidth=0.5)+
      annotate(geom="text",x=-Inf,y=Inf,
               hjust=-0.15,
               vjust=1+0.15*3,
               label=sprintf("λ = %.2f",LAMBDA),#sprintf("λ = %.2f",LAMBDA) #substitute(lambda==L,list(L=LAMBDA))
               parse=TRUE,
               size=5)+
      labs(x=expression("Expected -log"[10]~"(P-value)"),
           y=expression("Observed -log"[10]~"(P-value)"))+
      ggtitle(PLOTTITLE)+
      myggplottheme
    
  }else{
    
    #Format data.
    DATA<-data.frame(Observed=sort(PVALUES),
                     Expected=ppoints(N),
                     CI_lower=qbeta(p=(1-CONFINTERVAL)/2,shape1=1:N,shape2=N:1),
                     CI_upper=qbeta(p=(1+CONFINTERVAL)/2,shape1=1:N,shape2=N:1))
    if(MAXNPOINTS<Inf & nrow(DATA)>MAXNPOINTS){ #Sample points to be plotted.
      DATA<-DATA[sample(1:nrow(DATA),MAXNPOINTS,replace=FALSE),]
    }
    
    p<-ggplot(DATA)+
      geom_ribbon(mapping=aes(x=Expected,ymin=CI_lower,ymax=CI_upper),
                  alpha=0.1)+
      geom_segment(aes(x=0,y=SIGTHRESHOLD,xend=SIGTHRESHOLD,yend=SIGTHRESHOLD),
                   color="red",linetype=3,linewidth=0.5,alpha=0.5)+
      geom_segment(aes(x=SIGTHRESHOLD,y=0,xend=SIGTHRESHOLD,yend=SIGTHRESHOLD),
                   color="red",linetype=3,linewidth=0.5,alpha=0.5)+
      geom_point(aes(Expected,Observed),shape=16,size=1)+
      geom_abline(intercept=0,slope=1,color="blue",alpha=0.5)+
      geom_line(aes(Expected,CI_upper),linetype=3,linewidth=0.5)+
      geom_line(aes(Expected,CI_lower),linetype=3,linewidth=0.5)+
      annotate(geom="text",x=-Inf,y=Inf,
               hjust=-0.15,
               vjust=1+0.15*3,
               label=sprintf("λ = %.2f",LAMBDA),#sprintf("λ = %.2f",LAMBDA) #substitute(lambda==L,list(L=LAMBDA))
               parse=TRUE,
               size=5)+
      labs(x=expression("Expected P-value"),
           y=expression("Observed P-value"))+
      ggtitle(PLOTTITLE)+
      myggplottheme
    
  }
  
  return(p)
  
}

#PVALUES<-c(runif(1000,min=0,max=1),runif(500,min=0.9,max=1))
#PVALUES<-c(runif(1000,min=0,max=1),runif(500,min=0,max=0.01))
#plot_qqplot_pvalues(PVALUES=PVALUES, #Vector of p values.
#                    CONFINTERVAL=0.95,
#                    SIGTHRESHOLD=0.01,
#                    LOG10SCALE=TRUE,
#                    PLOTTITLE="QQ plot of p-values",
#                    MAXNPOINTS=Inf)
