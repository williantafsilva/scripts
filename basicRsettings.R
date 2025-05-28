############################################################################
################################ R SCRIPT ##################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

############################################################################
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

myggplottheme_nolegend<-theme(title=element_text(size=10,face="bold"),
                              axis.title=element_text(size=10,face="bold"),
                              axis.text=element_text(size=10),
                              axis.text.x=element_text(angle=60,size=8,vjust=0.5),
                              legend.position="none",
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

myggplottheme_blank_nolegend<-theme(title=element_text(size=10,face="bold"),
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
                                    panel.grid.major=element_blank(),
                                    panel.grid.minor=element_blank(),
                                    strip.background=element_rect(colour="black",
                                                                  fill="white"))
                           
############################################################################


















