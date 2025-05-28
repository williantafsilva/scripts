############################################################################
################################ R FUNCTION ################################
############################################################################
#Author: Willian T.A.F. Silva (willian.silva@evobiolab.com).
############################################################################

#Create a polygon with points, border and spikes. Radii must consist of positive values.
mypolygon<-function(radii=c(1,2,3,4,5),
                    color_points=c("red","red","red","red","red"),
                    color_border=c("blue","blue","blue","blue","blue"),
                    color_spikes=c("black","black","black","black","black"),
                    cex_points=0.5,lwd_border=1,lwd_spikes=1){
  if(length(color_points)!=length(radii)){color_points<-rep("red",length(radii))}
  if(length(color_border)!=length(radii)){color_border<-rep("blue",length(radii))}
  if(length(color_spikes)!=length(radii)){color_spikes<-rep("black",length(radii))}
  plot(0,type="n",axes=F,xlab="",ylab="",
       xlim=c(-max(radii),max(radii)),
       ylim=c(-max(radii),max(radii)))
  xvertices<-c()
  yvertices<-c()
  angles<-seq(from=0,to=2*pi,length.out=length(radii)+1)
  for(i in 1:length(radii)){
    xvertices<-c(xvertices,radii[i]*cos(angles[i]))
    yvertices<-c(yvertices,radii[i]*sin(angles[i]))
  }
  for(i in 1:length(xvertices)){
    lines(c(0,xvertices[i]),c(0,yvertices[i]),col=color_spikes[i],lwd=lwd_spikes)
    points(xvertices[i],yvertices[i],col=color_points[i],pch=16,cex=cex_points)
  }
  for(i in 1:length(xvertices)-1){
    lines(c(xvertices[i],xvertices[i+1]),c(yvertices[i],yvertices[i+1]),
          col=color_border[i],lwd=lwd_border)
    #polygon(x=c(xvertices[i],xvertices[i+1]),
    #        y=c(yvertices[i],yvertices[i+1]),
    #        border=color_border)
  }
  lines(c(xvertices[length(xvertices)],xvertices[1]),c(yvertices[length(xvertices)],yvertices[1]),
        col=color_border[length(xvertices)],lwd=lwd_border)
  #polygon(x=c(xvertices[length(xvertices)],xvertices[1]),
  #        y=c(yvertices[length(yvertices)],yvertices[1]),
  #        border=color_border)
}