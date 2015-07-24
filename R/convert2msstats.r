#' convert to MSstats
#' @export
#' @description convert msexp to msstats compatible output
#' @examples
#' data(feature_alignment_requant)
#' SDat = loadTransitonsMSExperiment(feature_alignment_requant)
#' msstatsformat = convert2MSstats(SDat2)
#' @param msexp to convert
convert2MSstats=function(msexp){
  #merge Intensity and pepinfo dataframes
  int = msexp$Intensity[match(rownames(msexp$Intensity), rownames(msexp$pepinfo)),]
  mergedDS = cbind(msexp$pepinfo, int)
  # remove spurious columns
  mergedDS = mergedDS[,!(colnames(mergedDS) %in% c("id","decoy"))]
  # melt the data into MSstats format
  meltedDS=melt(mergedDS,id=c("transition_group_id","aggr_Fragment_Annotation","PeptideSequence","PrecursorCharge","ProteinName"))
  meltedDS = meltedDS[,-match("transition_group_id",colnames(meltedDS))]
  colnames(meltedDS)[colnames(meltedDS) == "aggr_Fragment_Annotation"] = "FragmentIon"
  colnames(meltedDS)[colnames(meltedDS) == "value"] = "Intensity"
  colnames(meltedDS)[colnames(meltedDS) == "variable"] = "Run"
  meltedDS$"ProductCharge" = rep(1,dim(meltedDS)[1])
  meltedDS$"IsotopeLabelType" = rep("L",dim(meltedDS)[1])
  return(meltedDS)
}
