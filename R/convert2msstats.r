#' convert to MSstats
#' @export
#' @description convert msexp to msstats compatible output
#' @examples
#' NULL
convert2MSstats=function(msexp){
  #merge Intensity and pepinfo dataframes
  int = msexp$Intensity[match(rownames(msexp$Intensity), rownames(msexp$pepinfo)),]
  mergedDS = cbind(msexp$pepinfo, int)
  # remove spurious columns
  mergedDS = mergedDS[,!(colnames(mergedDS) %in% c("id","decoy"))]
  # melt the data into MSstats format
  meltedDS=melt(mergedDS,id=c("transition_group_id","aggr_Fragment_Annotation","PeptideSequence","PrecursorCharge","ProteinName"))
  colnames(meltedDS)[colnames(meltedDS) == "aggr_Fragment_Annotation"] = "FragmentIon"
  colnames(meltedDS)[colnames(meltedDS) == "value"] = "Intensity"
  colnames(meltedDS)[colnames(meltedDS) == "variable"] = "OrigFilename"
  return(meltedDS)
}
