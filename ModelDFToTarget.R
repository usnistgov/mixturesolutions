return(generalflattarget(m1,returndata=retdf,...)) #use the matrix to make a plot.  Return ggplot object
generalflattarget<-function(dfmatrix,flat=TRUE,componentnames=c("Brain","Liver","Placenta"),truemixvals=data.frame(Mix1=c(.25,.25,.5),Mix2=c(.5,.25,.25)),returndata=FALSE){
  require(dplyr)
  require(ggplot2)
  if(!"Mix"%in%colnames(dfmatrix)){
    dfmatrix<-as.data.frame(dfmatrix) #this input should* just be a data frame of the mixes and associated metadata
    dfmatrix[componentnames[1]]<-dfmatrix[componentnames[1]]-truemixvals[1,"Mix1"] #this could be looped i suppose?
    dfmatrix[componentnames[2]]<-dfmatrix[componentnames[2]]-truemixvals[2,"Mix1"]
    dfmatrix[componentnames[3]]<-dfmatrix[componentnames[3]]-truemixvals[3,"Mix1"]
    dfmatrix[,4]<-dfmatrix[,4]-truemixvals[1,"Mix2"]
    dfmatrix[,5]<-dfmatrix[,5]-truemixvals[2,"Mix2"]
    dfmatrix[,6]<-dfmatrix[,6]-truemixvals[3,"Mix2"]
    newnames<-colnames(dfmatrix)[4:6]#this is so not general
    #  tdfmatrix<-suppressMessages(melt(dfmatrix,measure.vars = c(componentnames,newnames),id.vars = colnames(dfmatrix)[colnames(dfmatrix)%in%c("Lab","Technology","Instrument","Round","Replicate","Plat")]))
    tdfmatrix<-suppressMessages(melt(dfmatrix,measure.vars = c(componentnames,newnames),id.vars = colnames(dfmatrix)[colnames(dfmatrix)%in%c("Lab","Technology","Instrument","Round","Replicate","Plat")]))
    require(dplyr)
    dfmatrix<-filter(tdfmatrix,!variable%in%newnames)
    dfmatrix$value2<-tdfmatrix[tdfmatrix$variable%in%newnames,"value"]
    dfmatrix$dfm<-sqrt(dfmatrix$value^2+dfmatrix$value2^2)
  }
  if(returndata==TRUE){return(dfmatrix)}
  g<-ggplot(data=dfmatrix)+
    geom_path(data=data.frame(circleFun(center=c(0,0),diameter=0.03,npoints=25),variable=componentnames[3]),aes(x,y),col="grey")+
    geom_path(data=data.frame(circleFun(center=c(0,0),diameter=0.06,npoints=25),variable=componentnames[2]),aes(x,y),col="grey")+
    geom_path(data=data.frame(circleFun(center=c(0,0),diameter=0.09,npoints=25),variable=componentnames[1]),aes(x,y),col="grey")+
    geom_point(data=data.frame(x=c(0,0,0),y=c(0,0,0),variable=componentnames),aes(x,y),size=3)+
    geom_point(data=dfmatrix,aes(x=value,y=value2,col=variable),size=5)+    scale_colour_manual(values=c("#6699FF","#99CC99","#CC6666"))+theme(aspect.ratio=1)+theme(legend.position="none")+xlab("")+ylab("")+theme(strip.background=element_rect(fill="white"))
  return(g)
}
#makes a flat (distance from target) target plot (ggplot object which can be modified)

