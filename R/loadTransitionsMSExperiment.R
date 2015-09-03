.syncrownames = function(msexp){
  for(fields in names(msexp)){
    tmp = msexp[[fields]]
    msexp[[fields]] = tmp[order(rownames(tmp)),]
  }
  return(msexp)
}
#' load msexperiment with nrt transtions and peptides
#' 
#' @description Selects top nrt transitions based on median transition intensity in all runs.
#' Selects top nr peptides based on median peptide intensity in all runs.
#' @export
#' @examples
#' library(imsbInfer)
#' data( feature_alignment_requant )
#' 
#' #SpecLib = ("C:/Users/witek/Google Drive/DataAnalysis/EBhardt/data/E1404301658-sample-SpecLib/feature_alignment_requant.tsv")
#' #SpecLib = ("/media/witold/data/googledrive/DataAnalysis/EBhardt/data/E1404301658-sample-SpecLib/feature_alignment_requant.tsv")
#' #obj  = fread(SpecLib)
#' #nrt = 3
#' #peptop = 3
#' obj =feature_alignment_requant
#' x = loadTransitonsMSExperiment( feature_alignment_requant , nrt= 3, peptop=3)
#' x = loadTransitonsMSExperiment( feature_alignment_requant )
#' table(table(feature_alignment_requant$transition_group_id))
#' table(table(x$pepinfo$transition_group_id))
#' x2 = loadTransitonsMSExperiment(feature_alignment_requant, nrt= 20, peptop=20)
#' table(table(x2$pepinfo$ProteinName))
#' table(table(x2$pepinfo$PeptideSequence))
#' dim(x)
#' head(x$pepinfo)
#' mypairs(x$Intensity[,1:3])
#' #check that ordering is consistent
#' xx = split2table(rownames(x$pepinfo),split="-")
#' stopifnot(xx[,1] == x$pepinfo$transition_group_id)
#' stopifnot(xx[,2] == x$pepinfo$aggr_Fragment_Annotation)
loadTransitonsMSExperiment = function(obj, nrt =3, peptop = 3){
  obj = prepareDF(obj)
  loccoll = Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "C")
  print(loccoll)
  ptm <- proc.time()
  cat("reading extended peptide information (creating msexperiment)\n - please be patient it make take a while (minutes)\n")
  
  msexp = read2msExperiment(obj)
  
  # long running function
  cat("extracting single transtion intensities\n - please be patient it make take a while (minutes)\n")
  transitiondata =  transitions2wide(obj)
  # this will read in also the full annotation (which peptide belongs to which protein)
  #rm(obj)
  gc()
  
  ##### selecting top 2-n fragments ####
  # long running
  cat("selecting top :", nrt , " transitions\n - please be patient it make take a while (minutes)\n")
  toptrans = selectTopFragmentsPerPeptide(transitiondata , nrt=nrt)
  gc()
  ##### 
  cat("aggregating peptide intensities based on top :", nrt , " transitons.\n")
  agrpeptide = .aggregatepeptide(toptrans)
  dim(agrpeptide)
  
  ## update the intensities with new intensities computed from top 2 transitions
  msexp$Intensity = agrpeptide[,2:dim(agrpeptide)[2]]
  rownames(msexp$Intensity) = agrpeptide$transition_group_id
  msexp$Intensity = msexp$Intensity[order(rownames(msexp$Intensity)),]
  
  msexp = .syncrownames(msexp)
  stopifnot(dim(msexp$Intensity) == dim(msexp$mz))
  stopifnot(rownames(msexp$Intensity) ==  rownames(msexp$mz))
  
  stopifnot(msexp$pepinfo$transition_group_id == rownames(msexp$Intensity))
  
  # select top peptides
  cat("selecting top :", peptop, " peptides per protein\n")
  toppep = selectTopPeptidesPerProtein(msexp ,peptop=peptop)
  stopifnot(rownames(toppep$Intensity) == rownames(toppep$pepinfo))
  
  # get the transitions belonging to the top peptides
  # select the toptransitions of the top peptides
  toptrans = toptrans[toppep$pepinfo$transition_group_id]
  msexp = toppep
  
  # create msexperiment containing transtions
  msExpTransition = function(toptrans,msexp){
    #select two columns only
    tt = toptrans[,transition_group_id,aggr_Fragment_Annotation]
    setkey(tt,transition_group_id,aggr_Fragment_Annotation)
    setkey(toptrans,transition_group_id,aggr_Fragment_Annotation)
    
    msexp$pepinfo = merge(as.data.frame(tt),msexp$pepinfo,by="transition_group_id")
    # R-FIX you need them as characters for rownames access
    msexp$pepinfo$transition_group_id = as.character(msexp$pepinfo$transition_group_id)
    
    stopifnot(msexp$pepinfo$transition_group_id==toptrans$transition_group_id)
    
    newkey = paste(msexp$pepinfo$transition_group_id,msexp$pepinfo$aggr_Fragment_Annotation,sep="-")
    rownames(msexp$pepinfo) = newkey
    
    msexp$rt = as.matrix(msexp$rt[msexp$pepinfo$transition_group_id , ])
    class(msexp$pepinfo$transition_group_id)
    msexp$rt = as.matrix(msexp$rt[msexp$pepinfo$transition_group_id , ])
    
    rownames(msexp$rt) = newkey
    msexp$score = as.matrix(msexp$score[msexp$pepinfo$transition_group_id , ])
    
    rownames(msexp$score) = newkey
    msexp$mz = as.matrix(msexp$mz[msexp$pepinfo$transition_group_id , ])
    rownames(msexp$mz) = newkey
    
    stopifnot(msexp$pepinfo$aggr_Fragment_Annotation==toptrans$aggr_Fragment_Annotation)    
    
    msexp$Intensity = as.matrix( toptrans[,3:dim(toptrans)[2],with=FALSE])
    rownames(msexp$Intensity) = newkey
    return(msexp)
  }
  
  msexp = msExpTransition(toptrans,toppep)
  cat("processing done in : ",proc.time() - ptm," [s] \n")
  gc()
  Sys.setlocale("LC_COLLATE", loccoll)
  print(Sys.getlocale("LC_COLLATE"))
  return(msexp)
}
