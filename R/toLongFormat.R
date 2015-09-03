#' extract transition intensities from input file. WARNING - slow running function
#' @description
#' the input table required fields are:
#' transition_group_id, 
#' aggr_Peak_Area,
#' aggr_Fragment_Annotion,
#' align_origfilename, 
#' @export
#' @examples
#' data(feature_alignment_requant)
#' far = feature_alignment_requant
#' df = prepareDF(feature_alignment_requant)
#' colnames(df)
#' tmp = transitions2wide(df)
#' names(df)
#' dim(feature_alignment_requant)
#' dim(tmp)[1]/dim(feature_alignment_requant)[1]*3
transitions2wide = function(far){
  loccoll = Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "C")
  
  #far = feature_alignment_requant
  ids = as.character(far$transition_group_id)
  apa = as.character(far$aggr_Peak_Area)
  afa = as.character(far$aggr_Fragment_Annotation)
  orig = basename(as.character(far$align_origfilename))
  orig = sub("_with_dscore.*","",orig)
  
  # split transition intensities
  transints = lapply(apa,function(x){unlist(strsplit(x,";"))})
  # split transition names
  transids = lapply(afa,function(x){unlist(strsplit(x,";"))})
  # prepare output
  lx = length(transids)
  transids <<-transids
  idss = vector(length=lx,mode="list")
  origf = vector(length=lx,mode="list")
  
  # extend
  for(i in 1:lx){
    l =  length(transids[[i]])
    idss[[i]] = rep(ids[i],l) # transition group ids
    origf[[i]] = rep(orig[i],l) # orig file name
  }
  # extend to wide format
  tmp = data.frame(transition_group_id = as.character(unlist(idss)),
                   align_origfilename = as.character(unlist(origf)),
                   aggr_Peak_Area = as.numeric(unlist(transints)),
                   aggr_Fragment_Annotation =  as.character(unlist( transids) ) )
  #got wide format here...
  data = dcast(tmp, transition_group_id + aggr_Fragment_Annotation ~ align_origfilename , value.var="aggr_Peak_Area")
  Sys.setlocale("LC_COLLATE", loccoll)
  
  return(data)
}