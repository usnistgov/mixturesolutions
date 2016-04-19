calcmrnafracgeneral <-
function(dat,spikeID="ERCC-",spikemassfraction=.1){
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
}
