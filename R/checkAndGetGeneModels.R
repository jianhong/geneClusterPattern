#' check and run getGeneModel
#' @noRd
#' @param genesList A list of GRanges.
#' @param ids IDs to be plotted. See \link{getGeneCluster}.
#' @param max_gap The maximal gaps from the first ID in `ids`.
#' @return a list with geneModel GRanges, region for that species,
#'  id0_width_pct and id_width_pct.
checkAndGetGeneModels <- function(genesList, ids, max_gap){
  stopifnot('genesList must have names'=
              length(names(genesList))==length(genesList))
  null <- lapply(genesList, function(.ele){
    if(!is(.ele, 'GRanges')){
      stop('genesList must be a list of GRanges', call. = FALSE)
    }
    if(length(.ele$gene_name)!=length(.ele)){
      stop('Elements in genesList must have metadata "gene_name"',
           call. = FALSE)
    }
    if(length(names(.ele))!=length(.ele)){
      stop('Elements in genesList must have names',
           call. = FALSE)
    }
  })
  
  return(lapply(genesList, getGeneModel,
                ids=ids, max_gap=max_gap))
}