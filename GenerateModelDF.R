preptargetplot<-function(mrnatype="consensus",splittype="none",modeltype="2mix",indf,retdf=FALSE,idvars=1,trueproportions=data.frame(mix1=c(.25,.25,.5),mix2=c(.5,.25,.25)),componentnames=c("Brain","Liver","Placenta"),mixnames=c("Mix1","Mix2"),...){
  indf<-normalizedf(indf,idvars) #apply a calcnormfactors (uqn) normalization
  #perhaps if we start splitting this into 5 separate data frames already it'd be better?
  #example:
  #C1<-indf[,grep(componentnames[1],colnames(indf))]
  #C2<-indf[,grep(componentnames[2],colnames(indf))]
  #C3<-indf[,grep(componentnames[3],colnames(indf))]
  #C4<-indf[,grep(componentnames[4],colnames(indf))]
  #C5<-indf[,grep(componentnames[5],colnames(indf))]
  #M1<-indf[,grep(mixnames[1],colnames(indf))]
  #M2<-indf[,grep(mixnames[2],colnames(indf))]
  #Some error checking here to make sure that the grep list is unique (no data is in more than 1 place) and complete (all mixes and components exist)

  normalizedf<-function(indf,idvars){
  require(edgeR)
  indf<-as.data.frame(indf) #in case it's actually read in as a data table.
  indf<-indf[rowSums(indf[-idvars])>0,]   #need to get rid of empty rows for upperquartile
  indf[-idvars]<-round(indf[-idvars]) #need to get rid of noninteger values for calcnormfactors
  indf[is.na(indf)]<-0  #this doesn't support NAs in the idvars.  hope that doesn't become an issue.
  normfac<-calcNormFactors(indf[,-idvars],method="upperquartile")
  normdf<-indf[-idvars] #subset to only the 'count' columns of df
  IT<-0 #initialize counter
  for(I in normfac){IT<-IT+1;normdf[,IT]<-normdf[,IT]*I} #calculate normalized counts by multiplying by normfac
  indf[-idvars]<-normdf #replace the counts with normalized counts
  return(indf)

}#applies calcnormfactors normalization: This looks good right now.
mfrac<-getmfrac(indf,mrnatype) #calculate the measured fraction ; do the
        getmfrac<-function(indf,type=mrnatype){
  #call appropriate helper functions based on the type of mRNA calculation i choose
  switch(type,
         internalconsensus = doubleconsensus(indf,modeltype),
         externalagreement = lookupmfrac(mixtureid),
         none=c(1,1,1),
         ercc=selectcomponents(calcmrnafractiongeneral(indf,"ERCC-")) #selectcomponents needs to be created. calcmrnafrac returns
         #each column's mfrac (need to select only the relevant columns from that data & normalize to 1)
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
minmodelsolve<-function(ssx=indf,mifrac,componentnames=c("Brain","Liver","Placenta"),mixnames=c("Mix1","Mix2")){
  modelproportions<-NULL
  #somewhat awkwardly, this assumes that the data consists only of 1 replicate of each component/mix.
  #i suppose i could do some pre-summarizing inside here?  Better: use a 'summarize' subroutine,elsewhere (in normalization?) to do that.


  for(idx in 1:length(mixnames)){
    form<-paste0(mixnames[idx],"~","I(",componentnames[1],"*",mifrac[1],")+","I(",componentnames[2],"*",mifrac[2],")+I(",componentnames[3],"*",mifrac[3],")+0")
    #the above formula is only going to work for 3component mixes, but it could be fixed for 2component
            modelproportions<-rbind(modelproportions,c(
      coefficients(lm(data=ssx,form))/
        sum(coefficients(lm(data=ssx,form))),
      mixnames[idx]))
  }
  rownames(modelproportions)<-c(componentnames,"mix")
  return(modelproportions)
} #test cases:  2mix_3compBLM: ok, but isnt returning the values of idvars how i would like.   2mix_2comp: Needs updating of the formula line.   2mix_3comp: untested.

magicmodelsolve<-function(ssx=indf,mifrac,componentnames=c("Mix1","Mix2","Mix3"),mixnames=c("Mix1","Mix2","Mix3"),idvars){
  labround<-NULL
  #awkwardly, this assumes that the data consists only of 1 replicate of each component/mix.

  for(idx in 1:length(mixnames)){
    labround<-rbind(labround,c(
      #Mm=c(signif(coefficients(lm(data=indf,I(mmix*valuemf[2])~I(lmmix*valuemf[3])+I(lmix*valuemf[1])+0))/sum(coefficients(lm(data=indf,I(mmix*valuemf[2])~I(lmmix*valuemf[3])+I(lmix*valuemf[1])+0))),digits=3)[1],0,signif(coefficients(lm(data=indf,I(mmix*valuemf[1])~I(lmmix*valuemf[3])+I(lmix*valuemf[2])+0))/sum(coefficients(lm(data=indf,I(mmix*valuemf[1])~I(lmmix*valuemf[3])+I(lmix*valuemf[2])+0))),digits=3)[2]),

      coefficients(lm(data=ssx,mixnames[idx]~I(componentnames[1]*mifrac[1])+I(componentnames[2]*mifrac[2])+I(componentnames[3]*mifrac[3])+0))/
        sum(coefficients(lm(data=ssx,mixnames[idx]~I(componentnames[1]*mifrac[1])+I(componentnames[2]*mifrac[2])+I(componentnames[3]*mifrac[3])+0))),
      mixnames[idx],idlist))
  } #I could stand to build this out for different numbers of components and mixes, but i'm not sure how.
  return(labround)
}#buildlater
magicv1modelsolve<-function(){
  #example:  L3=c(signif(coefficients(lm(data=indf,I(MixL3/mfrac["MixL3"]*valuemf[2])~lmmix+I(mmix/valuemf[2]*valuemf[1])+0))/sum(coefficients(lm(data=indf,I(MixL3/mfrac["MixL3"]*valuemf[2])~lmmix+I(mmix/valuemf[2]*valuemf[1])+0))),digits=3),0),

} #buildlater
calcmrnafracgeneral<-function(dat,spikeID="ERCC-",spikemassfraction=.1,componentnames){
  countcolumns<-which(unname(unlist(lapply(dat,class))=="numeric"))  #0) Identify which columns are counts and which are 'annotation':
  annotcolumn<-which(unname(unlist(lapply(dat,class))!="numeric"))
  #1) Identify which counts are Spike-In and which are not
  ercc<-rownames(dat)[(grep(spikeID,rownames(dat)))]   #one way to identify spikes, if row names = spikeID
  if(length(ercc)==0){ercc<-grep(spikeID,dat[,annotcolumn[1]])} #Try a different way, assuming that the name is in the first annotation column...
  if(length(ercc)==0){stop("I can't identify the spike-ins within your count table.   The spikeID variable should be set to something which uniquely identifies spike-ins.   The SpikeID should be found in either Rownames or the FIRST non-numeric column.")}
  nonercc<-!(1:length(dat[,countcolumns[1]]))%in%ercc
  count<-rbind(colSums(dat[ercc,countcolumns]),colSums(dat[nonercc,countcolumns])) #determines the counts for spikes and non-spikes.
  ercc.targ<-spikemassfraction  #defines the "targeted" mass fraction for spikes : Either a vector with length = #columns,or a scalar
  mRNA.frac<-ercc.targ*count[2,]/count[1,]  #calculates an mRNA fraction based on those available data
  #return(mRNA.frac)  #this part doesn't normalize to one, but that's not exactly complicated.
  return(mRNA.frac[names(mRNA.frac)%in%componentnames]/sum(mRNA.frac[names(mRNA.frac)%in%componentnames])) #working on BLM.
}
#test cases:  2mix_3compBLM (working), 2mix_3comp (NA,nospikes), 2mix_2comp (untested, unsummarized), 3mix(build later)



doubleconsensus<-function(dataf,modeltype){
  fitmodel<-function(indf,modeltype){
    switch(modeltype,
           twomix=minmodelsolve(indf,mifrac,componentnames,mixnames,idvars),
           threemix=magicmodelsolve(indf,mifrac,componentnames,mixnames,idvars))
    return(data.frame(coef1=coefficients(fit)/sum(coefficients(fit)),coef2=coefficients(fit2)/sum(coefficients(fit2))))
  }
  f2opt<-function(par,indata){sum(abs(fitmodel(mifrac=c(par,1),indata)[,1]-trueproportions[,1]),abs(fitmodel(mifrac=c(par,1),indata)[,2]-trueproportions[,2]))}
  opted<-optim(par=c(1,1),f2opt,indata=dataf)$par
  return(c(opted,1)/sum(c(opted,1)))
}#test cases:  2mix_3compBLM (untested), 2mix_3comp (untested), 2mix_2comp (untested, unsummarized), 3mix(build later)
