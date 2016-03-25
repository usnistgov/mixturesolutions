preptargetplot<-function(mrnatype="internalconsensus",splittype="none",prenormalized=TRUE,modeltype="twomix",indf,retdf=FALSE,idvars=1,trueproportions=data.frame(mix1=c(.25,.25,.5),mix2=c(.5,.25,.25)),componentnames=c("Brain","Liver","Placenta"),mixnames=c("Mix1","Mix2"),...){
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

  }

  if(prenormalized==FALSE){indf<-normalizedf(indf,idvars)} #apply a calcnormfactors (uqn) normalization
  ###perhaps if we start splitting this into 5 separate data frames already it'd be better?
  #example of above:
  #C1<-indf[,grep(componentnames[1],colnames(indf))]
  #C2<-indf[,grep(componentnames[2],colnames(indf))]
  #C3<-indf[,grep(componentnames[3],colnames(indf))]
  #C4<-indf[,grep(componentnames[4],colnames(indf))]
  #C5<-indf[,grep(componentnames[5],colnames(indf))]
  #M1<-indf[,grep(mixnames[1],colnames(indf))]
  #M2<-indf[,grep(mixnames[2],colnames(indf))]

  ###I would also like to make a pass at turning all of the arguments to this function into a single object so they can get passed down more easily

  ###Some error checking here to make sure that the grep list is unique (no data is in more than 1 place) and complete (all mixes and components exist)

 #applies calcnormfactors normalization: This looks good right now.
  getmfrac<-function(indf,type=mrnatype){
    #call appropriate helper functions based on the type of mRNA calculation i choose
    mfrac<-switch(type,
           internalconsensus = doubleconsensus(indf,modeltype,mixnames,componentnames,trueproportions),
           externalagreement = lookupmfrac(mixtureid),
           none=c(1,1,1),
           ercc=calcmrnafracgeneral(indf,"ERCC-"),componentnames = componentnames[selectcomponents(calcmrnafracgeneral(indf,"ERCC-"),componentnames = componentnames)] #selectcomponents needs to be created. calcmrnafrac returns
           #each column's mfrac (need to select only the relevant columns from that data & normalize to 1)
    )
    mfrac<-mfrac/sum(mfrac)
    return(mfrac)
  }#call mfraction-calculating helper functions based on the type of measuredRNA
  mfrac<-getmfrac(indf,mrnatype) #calculate the measured fraction

        #calculation i choose
  splitmultidf<-function(indf,splittype,splitcolumn,idvars){
    switch(splittype,
           none=return(indf),
           globalcon = makelabrounda() #i need to split the 'makelabround' into a helper function that splits and a secondary
           #function that solves the model
           #the helper should just return a list of the dataframes & a list of the mfracs?
           #then the model-solving function can 'just' determine if it's a singledf or multidf
           #and return output accordingly.
    )
  } #converts a dataframe with multiple count tables into individuals
  #splitting types need work ; default to 'none' atm.  -- perhaps i need to stop even trying to split from bigdfs
  splitm1<-splitmultidf(indf,splittype,splitcolumn) #separate out multiple data files (optional)

  getm1<-function(indf,mfrac,type=modeltype){
    switch(type,
           twomix=minmodelsolve(ssx=indf,mfrac=mfrac,componentnames=componentnames,mixnames=mixnames),
           threemix=magicmodelsolve()

    )
  }
  m1<-getm1(indf,mfrac,modeltype) #use the model to create a matrix of deviation from target
#calls a model-solving helper function (eg: minmodelsolve)
return(m1)
}
#miscellaneous
minmodelsolve<-function(ssx=indf,mfrac,componentnames=c("Brain","Liver","Placenta"),mixnames=c("Mix1","Mix2")){
  modelproportions<-NULL
  #somewhat awkwardly, this assumes that the data consists only of 1 replicate of each component/mix.
  #i suppose i could do some pre-summarizing inside here?  Better: use a 'summarize' subroutine,elsewhere (in normalization?) to do that.

  for(idx in 1:length(mixnames)){
    #more generalized...
    form<-paste0(mixnames[idx],"~")
    for(cindex in (1:(length(componentnames)))){
       form<-paste0(form,"I(",componentnames[cindex],"*",mfrac[cindex],")+")
      }
    form<-paste0(form,0)

    modelproportions<-rbind(modelproportions,c(
      coefficients(lm(data=ssx,form))/
        sum(coefficients(lm(data=ssx,form))),
      mixnames[idx]))
  }
  colnames(modelproportions)<-c(componentnames,"mix")
  modelproportions<-as.data.frame(modelproportions) #converts proportions to numeric rather than character.
  for(I in 1:length(componentnames)){modelproportions[,I]<-as.numeric(as.character(modelproportions[,I]))}
  return(modelproportions)
} #test cases:  2mix_3compBLM: ok   2mix_2comp: untested, possible working 2mix_3comp: untested, probably fine
magicmodelsolve<-function(ssx=indf,mfrac,componentnames=c("Mix1","Mix2","Mix3"),mixnames=c("Mix1","Mix2","Mix3"),idvars){
  labround<-NULL
  #awkwardly, this assumes that the data consists only of 1 replicate of each component/mix.

  for(idx in 1:length(mixnames)){
    labround<-rbind(labround,c(
      #Mm=c(signif(coefficients(lm(data=indf,I(mmix*valuemf[2])~I(lmmix*valuemf[3])+I(lmix*valuemf[1])+0))/sum(coefficients(lm(data=indf,I(mmix*valuemf[2])~I(lmmix*valuemf[3])+I(lmix*valuemf[1])+0))),digits=3)[1],0,signif(coefficients(lm(data=indf,I(mmix*valuemf[1])~I(lmmix*valuemf[3])+I(lmix*valuemf[2])+0))/sum(coefficients(lm(data=indf,I(mmix*valuemf[1])~I(lmmix*valuemf[3])+I(lmix*valuemf[2])+0))),digits=3)[2]),

      coefficients(lm(data=ssx,mixnames[idx]~I(componentnames[1]*mfrac[1])+I(componentnames[2]*mfrac[2])+I(componentnames[3]*mfrac[3])+0))/
        sum(coefficients(lm(data=ssx,mixnames[idx]~I(componentnames[1]*mfrac[1])+I(componentnames[2]*mfrac[2])+I(componentnames[3]*mfrac[3])+0))),
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
}#test cases:  2mix_3compBLM (working), 2mix_3comp (NA,nospikes), 2mix_2comp (untested, unsummarized), 3mix(build later)


doubleconsensus<-function(dataf,modeltype,mixnames,componentnames,trueproportions){
  fitmodel<-function(indf,mfrac=c(1,1,1),modeltype,mixnames,componentnames){
    fit<-switch(modeltype,
           twomix=minmodelsolve(indf,mfrac,componentnames,mixnames),
           threemix=magicmodelsolve(indf,mfrac,componentnames,mixnames))
    return(fit)
  }
  f2opt<-function(par,indata,modeltype,mixnames,componentnames,trueproportions){sum(abs(fitmodel(mfrac=c(par,1),indata,modeltype=modeltype,mixnames=mixnames,componentnames=componentnames)[1,1:length(trueproportions[,1])]-trueproportions[,1]),
                                                     abs(fitmodel(mfrac=c(par,1),indata,modeltype=modeltype,mixnames=mixnames,componentnames=componentnames)[2,1:length(trueproportions[,2])]-trueproportions[,2]))}
  opted<-optim(par=c(1,1),f2opt,indata=dataf,modeltype=modeltype,mixnames=mixnames,componentnames=componentnames,trueproportions=trueproportions)$par
  return(c(opted,1)/sum(c(opted,1)))
}#test cases:  2mix_3compBLM (untested), 2mix_3comp (untested), 2mix_2comp (untested, unsummarized), 3mix(build later: this CLEARLY depends on trueproportions having only 2 columns)
selectcomponents<-function(dataf,componentnames){which(names(dataf)%in%componentnames)}

calcmrnafracgeneral<-function(dat,spikeID="ERCC-",spikemassfraction=.1){
  #1) Identify which counts are Spike-In and which are not
  #0) Identify which columns are counts and which are 'annotation':
  countcolumns<-which(unname(unlist(lapply(dat,class))=="numeric"))
  annotcolumn<-which(unname(unlist(lapply(dat,class))!="numeric"))
  ercc<-rownames(dat)[which(substr(rownames(dat),1,5)==spikeID)]           #one way to identify spikes, if row names = spikeID
  if(length(ercc)==0){ercc<-grep(spikeID,dat[,annotcolumn[1]])} #assuming that the name is in the first annotation column...
  if(length(ercc)==0){stop("I can't identify the spike-ins within your count table.   The spikeID variable should be set to something which uniquely identifies spike-ins.   Rownames are first checked for names, then if there are non-numeric columns, only the FIRST is checked for gene names. ")}
  nonercc<-!(1:length(dat[,countcolumns[1]]))%in%ercc

  count<-rbind(colSums(dat[ercc,countcolumns]),colSums(dat[nonercc,countcolumns])) #determines the counts for spikes and non-spikes.
  ercc.targ<-spikemassfraction  #defines the "targeted" mass fraction for spikes : Either a vector with length = #columns,or a scalar
  mRNA.frac<-ercc.targ*count[2,]/count[1,]  #calculates an mRNA fraction based on those available data
  #this part doesn't normalize to one, but that's not exactly complicated.
  return(mRNA.frac)
}#2mix_3comp_blm tested.