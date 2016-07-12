#' Calculate the fraction of mixture components in a mixture dataset.
#'
#' \code{mixmetrics} solves a linear model based on the components
#' and proportions listed.  Returns a modified object with model results and
#' predicted count values.
#'
#'
#' @param mrnatype Internalconsensus,externalagreement,none,or ercc;
#' How to calculate the measured fraction in a sample.
#' @param splittype 'none' - not currently implemented:
#' How to split the data rows based on columns in idvars.
#' @param prenormalized TRUE/FALSE - if data is already normalized.
#' @param modeltype "twomix","threemix":
#' How many different mixtures are present in your dataset.
#' @param indf Input dataframe:  counts and idvariables.
#' @param idvars List of columns (or names?) that are not data, but identifiers
#' of gene names(1st) or splitting features (not 1st, not implemented).
#' @param trueproportions Data frame consisting of the true mixture
#' proportions of each mix, in order.
#' @param componentnames List of column names of indf for each component.
#' @param idnames What you want to call your identifier field. (default:id)
#' @param mixnames List of column names of indf corresponding to each mix.
#' @param annot How to refer to the various copies of each component/mix
#' after they are split/merged.
#' @param ... Anything else.
#' @export mixmetrics
#'
#'
# I would like to update this to : (done): a) Return an object that's just a data structure with various bits
# b) assign a method to coef to extract the lm-fit coefficients to the output data structure
# c) assign a method to .predict to predict the counts from the output data structure
# d) Determine a way to piece together multiple datasets for iterative studies:
#    [currently: Just rbind multiple results]
mixmetrics <-function(mrnatype="internalconsensus",splittype="none",prenormalized=TRUE,modeltype="twomix",indf,idnames="id",idvars=1,trueproportions,componentnames=c("Brain","Liver","Placenta"),mixnames=c("Mix1","Mix2"),annot=c("Replicate"),...){
##internally-used function definition
  getmfrac<-function(indf,type=preps$mrnatype){
    #call appropriate helper functions based on the type of mRNA calculation i choose
    mfrac<-switch(type,
                  internalconsensus = doubleconsensus(indf,preps),#modeltype,mixnames,componentnames,trueproportions),
                  #externalagreement = lookupmfrac(mixtureid),
                  none=c(1,1,1),
                  ercc=calcmrnafracgeneral(indf[c(1,selectcomponents(indf,preps))])#
    )
    mfrac<-mfrac/sum(mfrac)
    return(mfrac)
  }#call mfraction-calculating helper functions based on the type of measuredRNA
  collapsereps<-function(indf,preps,annot){
    reslist<-NULL
    for(I in 1:length(eval(preps$componentnames))){
      if(length(grep(eval(preps$componentnames)[I],colnames(indf)))==0){warning(call. = TRUE,... = paste0("One of your componentnames, ",eval(preps$componentnames)[I]," couldn't be found in the data:  Do colnames match componentnames?"))}
      reslist<-c(reslist,grep(eval(preps$componentnames)[I],colnames(indf)))
    }
    for(I in 1:length(eval(preps$mixnames))){
      if(length(grep(eval(preps$mixnames)[I],colnames(indf)))==0){warning(call. = TRUE,... = paste0("One of your mixnames, ",eval(preps$mixnames)[I]," couldn't be found in the data:  Do colnames match mixnames?"))}
      reslist<-c(reslist,grep(eval(preps$mixnames)[I],colnames(indf)))}
    if(length(unique(reslist))!=length(reslist)){warning("One or more entities appear to have non-unique names.  Can you make sure all of your components and mixes have unique names?")}
    ###Some error checking here to make sure that the grep list is unique (no data is in more than 1 place) and complete (all mixes and components exist)
    C1<-indf[,c(grep(eval(preps$componentnames)[1],colnames(indf)))]; C1l<-NULL;for(I in 1:ncol(C1)){C1l<-rbind(C1l,data.frame(id=indf[,idvars],data=C1[,I],annot=I))};colnames(C1l)<-c(idnames,"data",annot);C1l$Source<-eval(preps$componentnames)[1]
    if(length(eval(preps$componentnames))>1){C2<-indf[,c(grep(eval(preps$componentnames)[2],colnames(indf)))];C2l<-NULL;for(I in 1:ncol(C2)){C2l<-rbind(C2l,data.frame(id=indf[,idvars],data=C2[,I],annot=I))};colnames(C2l)<-c(idnames,"data",annot);C2l$Source<-eval(preps$componentnames)[2]}
    if(length(eval(preps$componentnames))>2){C3<-indf[,c(grep(eval(preps$componentnames)[3],colnames(indf)))] ;C3l<-NULL;for(I in 1:ncol(C3)){C3l<-rbind(C3l,data.frame(id=indf[,idvars],data=C3[,I],annot=I))};colnames(C3l)<-c(idnames,"data",annot);C3l$Source<-eval(preps$componentnames)[3]}
    if(length(eval(preps$componentnames))>3){warning("You seem to have input more than 3 component names.  This behavior is not currently supported.")}
    if(length(eval(preps$componentnames))>3){C4<-indf[,c(grep(eval(preps$componentnames)[4],colnames(indf)))] ;C4l<-NULL;for(I in 1:ncol(C4)){C4l<-rbind(C4l,data.frame(id=indf[,idvars],data=C4[,I],annot=I))};colnames(C4l)<-c(idnames,"data",annot);C4l$Source<-eval(preps$componentnames)[4]}
    if(length(eval(preps$componentnames))>4){C5<-indf[,c(grep(eval(preps$componentnames)[5],colnames(indf)))] ;C5l<-NULL;for(I in 1:ncol(C5)){C5l<-rbind(C5l,data.frame(id=indf[,idvars],data=C5[,I],annot=I))};colnames(C5l)<-c(idnames,"data",annot);C5l$Source<-eval(preps$componentnames)[5]}
    M1<-indf[,c(grep(eval(preps$mixnames)[1],colnames(indf)))];M1l<-NULL;for(I in 1:ncol(M1)){M1l<-rbind(M1l,data.frame(id=indf[,idvars],data=M1[,I],annot=I))};colnames(M1l)<-c(idnames,"data",annot);M1l$Source<-eval(preps$mixnames)[1]
    if(length(eval(preps$mixnames))>1){M2<-indf[,c(grep(eval(preps$mixnames)[2],colnames(indf)))];M2l<-NULL;for(I in 1:ncol(M2)){M2l<-rbind(M2l,data.frame(id=indf[,idvars],data=M2[,I],annot=I))};colnames(M2l)<-c(idnames,"data",annot);M2l$Source<-eval(preps$mixnames)[2]}
    if(length(eval(preps$mixnames))>2){warning("You have input more than 2 mix names.  This behavior is not currently supported")}

    if(length(eval(preps$mixnames))>2){M3<-indf[,c(grep(eval(preps$mixnames)[3],colnames(indf)))];M3l<-NULL;for(I in 1:ncol(M3)){M3l<-rbind(M3l,data.frame(id=indf[,idvars],data=M3[,I],annot=I))};colnames(M3l)<-c(idnames,"data",annot);M3l$Source<-eval(preps$mixnames)[3]}
    if(length(eval(preps$mixnames))>3){M4<-indf[,c(grep(eval(preps$mixnames)[4],colnames(indf)))];M4l<-NULL;for(I in 1:ncol(M4)){M4l<-rbind(M4l,data.frame(id=indf[,idvars],data=M4[,I],annot=I))};colnames(M4l)<-c(idnames,"data",annot);M4l$Source<-eval(preps$mixnames)[4]}

    ###I need to do something to get the ID vector put as part of the dataset ;
    #add everything that *did* get created in the set of C1:C5&M1:M4 into the output list.
    #apparently the environment is changing pretty quickly around here so i need to define things before searching...If any additional changes occur below, they need to be pre-defined here!
    listout<-NULL;outdf<-NULL;IT<-0;allnames<-c(eval(preps$componentnames),eval(preps$mixnames));J<-0

    listout<-c(grep("^C[0-9]l$",ls()),grep("^M[0-9]l$",ls())) #find all of the objects we created fitting the pattern.
    #This pattern should be sufficiently rigid that no false positives show, but searching only within environmentwould help.  If i knew how to do that.
    if(length(listout)!=length(allnames)){stop("Something went horribly wrong in replicate collapsing.")}
    for(J in listout){
      IT<-IT+1
      outdf<-rbind(outdf,get(ls()[J]))} #paste them all together into an object.The type of it is the tricky part...
    #I chose to make it a (molten) data frame.  It still needs to be dcast for lm, though, and probably other things.
    outdf<-reshape2::dcast(outdf,get(idnames)~Source,fun.aggregate=mean,na.rm=TRUE,value.var = 'data') #This is going to have to work for now. It has issues.
    return(outdf)
  }
  normalizedf<-function(indf,idvars){
    indf<-as.data.frame(indf) #in case it's read in as a data table.
    indf[is.na(indf)]<-0  #this doesn't support NAs in the idvars.  hope that doesn't become an issue.
    indf<-indf[rowSums(indf[-idvars])>0,]   #get rid of empty rows for upperquartile
    indf[-idvars]<-round(indf[-idvars]) #get rid of noninteger values for calcnormfactors
    normfac<-edgeR::calcNormFactors(indf[,-idvars],method="upperquartile")
    normdf<-indf[-idvars] #subset to only the 'count' columns of df
    IT<-0 #initialize counter
    for(I in normfac){IT<-IT+1;normdf[,IT]<-normdf[,IT]*I} #calculate normalized counts by multiplying by normfac
    indf[-idvars]<-normdf #replace the counts with normalized counts
    return(indf)

  }  #applies calcnormfactors normalization: This looks good right now.
  cellmixmodelsolve<-function(indf,preps){
    #yeah you're going to lose all of the 0s. It also sucks for BLM.
    #I can't imagine a different way.
    # this can run as a 'blind solve N components' or as a 'known signatures for components'
    # but they both suck.
    indfb<-indf[rowSums(indf[preps$mixnames])>0,]
    indfb<-indfb[rowSums(indfb[preps$componentnames])>0,]
    return(coefficients(CellMix::ged(object=as.matrix(indfb[preps$mixnames]),x=length(preps$componentnames))))
  }
  splitmultidf<-function(indf,splittype,splitcolumn,idvars){
    switch(splittype,
           none=return(indf)#,
           #globalcon = makelabrounda() #i need to split the 'makelabround' into a helper function that splits and a secondary
           #function that solves the model
           #the helper should just return a list of the dataframes & a list of the mfracs?
           #then the model-solving function can 'just' determine if it's a singledf or multidf
           #and return output accordingly.
    )
  } #INCOMPLETE.  Converts a dataframe with multiple count tables into individuals
##end function definitions


  #Turning all of the arguments to this function into a single structure so they can get passed into sub-functions more easily
  #preps<-list(modeltype=modeltype,mrnatype=mrnatype,trueproportions=trueproportions,mixnames=mixnames,componentnames=componentnames,idvars=idvars,idnames=idnames)
  preps<-as.list(match.call())
  if(prenormalized==FALSE){indf<-normalizedf(indf,idvars)} #apply a calcnormfactors (uqn) normalization if needed
  if(length(grep(componentnames[1],colnames(indf)))>1){
    indf=collapsereps(indf,preps,annot)
  } #converts replicate (underscore-separated) data into annotated(with whatever was in the underscores...) data

  mfrac<-getmfrac(indf,preps$mrnatype) #calculate the measured fraction
  preps<-c(preps,list(mfrac=mfrac))#adds the mfrac to the data structure
  #splitm1<-splitmultidf(indf,splittype,splitcolumn) #separate out multiple data files (future)

  getm1<-function(indf,preps){
    switch(preps$modeltype,
           twomix=minmodelsolve(ssx=indf,preps=preps),
           #threemix=magicmodelsolve(),
           #would be nice to use cellmix as an option
           cellmix=cellmixmodelsolve(indf,preps)

    )
  }
  m1<-getm1(indf,preps=preps) #use the model to create a matrix of deviation from target
#calls a model-solving helper function (eg: minmodelsolve)
  preps<-c(preps,list(results=m1,indf=indf))
  targetdf<-data.frame(x=c(eval(preps$trueproportions)[,1],as.numeric(preps$results[1,1:length(eval(preps$componentnames))])),
                       y=c(eval(preps$trueproportions)[,2],as.numeric(preps$results[2,1:length(eval(preps$componentnames))])),
                       col=c(eval(preps$componentnames),eval(preps$componentnames)),
                       meta=c(rep("expected",length(eval(preps$trueproportions)[,1])),rep("observed",length(eval(preps$componentnames))))
                      )#builds a data frame that's standardized for converting into a target plot
  preps<-c(preps,list(targetdf=targetdf))
    return(preps)
}
