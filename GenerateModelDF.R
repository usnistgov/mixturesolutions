preptargetplot<-function(mirnatype="consensus",splittype="none",modeltype="2mix",indf,retdf=FALSE,idvars=1,componentnames=c("Brain","Liver","Placenta"),...){
indf<-normalizedf(indf,idvars) #apply a calcnormfactors (uqn) normalization

      normalizedf<-function(indf,idvars){
  require(edgeR)
  indf<-as.data.frame(indf) #in case it's actually read in as a data table.
  indf<-indf[rowSums(indf[-idvars])>0,]   #need to get rid of empty rows
  indf[-idvars]<-round(indf[-idvars])
  indf[is.na(indf)]<-0  #this doesn't support NAs in the idvars.  hope that doesn't become an issue.
  normfac<-calcNormFactors(indf[,-idvars],method="upperquartile") #UQN hates miRNA data because there're too many 0s.
  IT<-0
  normdf<-indf[-idvars]
  for(I in normfac){IT<-IT+1;normdf[,IT]<-normdf[,IT]*I}
  indf[-idvars]<-normdf
  return(indf)

}#applies calcnormfactors normalization
mifrac<-getmifrac(indf,mirnatype) #calculate the measured fraction ;
        getmifrac<-function(indf,type=mirnatype){
  #call appropriate helper functions based on the type of miRNA calculation i choose
  switch(type,
         consensus = doubleconsensus(indf),
         none=c(1,1,1),
         ercc=selectcomponents(calcmrnafractiongeneral(indf,"ERCC-")), #selectcomponents needs to be created. calcmrnafrac returns
         #each column's mfrac (needs to )
         #      labcon=,
         labround=alllabroundmfrac(indf),
         #       platform=
         roundrep = allreplicate(indf)
  )
}#call mfraction-calculating helper functions based on the type of measuredRNA
        #calculation i choose
splitm1<-splitmultidf(indf,splittype,splitcolumn) #separate out multiple data files (optional)
          splitmultidf<-function(indf,splittype,splitcolumn,idvars){
  switch(type,
         none=return(indf),
         globalcon = makelabrounda() #i need to split the 'makelabround' into a helper function that splits and a secondary
         #function that solves the model
         #the helper should just return a list of the dataframes & a list of the mfracs?
         #then the model-solving function can 'just' determine if it's a singledf or multidf
         #and return output accordingly.
  )
} #converts a dataframe with multiple count tables into individuals
          #splitting types need work ; default to 'none' atm.  -- perhaps i need to stop even trying to split from bigdfs
m1<-getm1(indf,mifrac,mirnatype) #use the model to create a matrix of deviation from target
getm1<-function(indf,mifrac,type=modeltype){
  switch(type,
         twomix=minmodelsolve(indf,mifrac),
         threemix=magicmodelsolve()

  )
}#calls a model-solving helper function (eg: minmodelsolve)
}
#miscellaneous
minmodelsolve<-function(ssx=indf,mifrac,componentnames=c("Brain","Liver","Placenta"),mixnames=c("Mix1","Mix2"),idvars){
  labround<-NULL
  #somewhat awkwardly, this assumes that the data consists only of 1 replicate of each component/mix.
  #i suppose i could do some pre-summarizing inside here?  Better: use the 'split' subroutine to do that.

  for(idx in 1:length(mixnames)){
    labround<-rbind(labround,c(
      coefficients(lm(data=ssx,mixnames[idx]~I(componentnames[1]*mifrac[1])+I(componentnames[2]*mifrac[2])+I(componentnames[3]*mifrac[3])+0))/
        sum(coefficients(lm(data=ssx,mixnames[idx]~I(componentnames[1]*mifrac[1])+I(componentnames[2]*mifrac[2])+I(componentnames[3]*mifrac[3])+0))),
      mixnames[idx],idlist))
  } #I could stand to build this out for different numbers of components and mixes, but i'm not sure how.
  return(labround)
} #untested
circleFun<- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  } #makes a circle for ggplot.
magicmodelsolve<-function(){

}
magicv1modelsolve<-function(){

}
makelabround<-function(mifrac=c(1,1,1),temp){
  labround<-NULL
  for(I in levels(as.factor(temp$Lab))){
    for(J in levels(as.factor(temp$Round))){
      for(K in levels(as.factor(temp$Instrument))){
        ssx<-subset(temp,Lab==I&Round==J&Instrument==K)
        if(length(ssx[,1])>0){
          if(length(ssx$variance)==0){
          labround<-rbind(labround,c(coefficients(lm(data=ssx,Mix1~I(Brain*mifrac[1])+I(Liver*mifrac[2])+I(Placenta*mifrac[3])+0))/sum(coefficients(lm(data=ssx,Mix1~I(Brain*mifrac[1])+I(Liver*mifrac[2])+I(Placenta*mifrac[3])+0))),"Mix1",I,J,K))
          labround<-rbind(labround,c(coefficients(lm(data=ssx,Mix2~I(Brain*mifrac[1])+I(Liver*mifrac[2])+I(Placenta*mifrac[3])+0))/sum(coefficients(lm(data=ssx,Mix2~I(Brain*mifrac[1])+I(Liver*mifrac[2])+I(Placenta*mifrac[3])+0))),"Mix2",I,J,K))
          }
          else{
            labround<-rbind(labround,c(coefficients(lm(data=ssx,Mix1~I(Brain*mifrac[1])+I(Liver*mifrac[2])+I(Placenta*mifrac[3])+0,weights=1/variance))/sum(coefficients(lm(data=ssx,Mix1~I(Brain*mifrac[1])+I(Liver*mifrac[2])+I(Placenta*mifrac[3])+0,weights=1/variance))),"Mix1",I,J,K))
            labround<-rbind(labround,c(coefficients(lm(data=ssx,Mix2~I(Brain*mifrac[1])+I(Liver*mifrac[2])+I(Placenta*mifrac[3])+0,weights=1/variance))/sum(coefficients(lm(data=ssx,Mix2~I(Brain*mifrac[1])+I(Liver*mifrac[2])+I(Placenta*mifrac[3])+0,weights=1/variance))),"Mix2",I,J,K)
          }
          }
      }}}
  labround<-as.data.frame(labround);  colnames(labround)<-c("Brain","Liver","Placenta","Mix","Lab","Round","Plat")
  for(I in 1:3){labround[,I]<-as.double(as.character(labround[,I]))}
  m1<-cbind(subset(labround,Mix=="Mix1"),subset(labround,Mix=="Mix2"))[,c(1:3,8:10,12,13,14)]
  colnames(m1)[4:6]<-c("B2","L2","P2")

    return(m1)
}