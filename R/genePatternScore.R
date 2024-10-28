#' Calculate the gene pattern score
#' @description
#' The gene pattern score is the normalized standard divination of 
#' edit (Levenshtein) distance of given ids from different species.
#' The maximal gene pattern score is 100.
#' Gene pattern score = 100 / ((sd(adist)+1)*k/length(ids))
#' @param genesList A list of GRanges.
#' @param ids IDs to be plotted. See \link{getGeneCluster}.
#' @param k The maximum number of nearest genes used in \link{getGeneCluster}.
#' @param max_gap The maximal gaps from the first ID in `ids`.
#' @param output Output values
#' @return the standard deviation or edit distance values.
#' @importFrom stats sd
#' @export
#' @examples
#' # example code
#' 
genePatternScore <- function(genesList, ids, k=length(ids), max_gap=1e7,
                             output=c("score", "sd", "distance")){
  output <- match.arg(output)
  stopifnot(is.numeric(k))
  stopifnot(length(k)==1)
  fac <- k/length(ids)
  geneIdMap <- c(seq(65, 90),# A-Z
                 seq(97, 122), # a-z
                 seq(48, 57), # 0-9
                 seq(32, 47), 
                 seq(58, 64),
                 seq(91, 96),
                 seq(123, 125))
  geneIdMap <- intToUtf8(geneIdMap)
  geneIdMap <- strsplit(geneIdMap, '')[[1]]
  
  if(length(ids)>length(geneIdMap)){
    stop('Too much ids. Maximal ', length(geneIdMap), ' ids are supported!')
  }
  
  # get the strings of gene pattern for each species
  geneModels <- checkAndGetGeneModels(genesList, ids=ids, max_gap=max_gap)
  geneIds <- lapply(geneModels, function(.ele) names(.ele$geneModel))
  
  # id map
  geneIdMap <- geneIdMap[seq_along(ids)]
  names(geneIdMap) <- ids
  
  ## mask the non-essential IDs
  maskedGeneIds <- lapply(geneIds, function(.ele){
    ch <- geneIdMap[.ele]
    ch[is.na(ch)] <- '~' # 126
    sub('^~+', '',
        sub('~+$', '',
            paste(ch, collapse=''))) # remove ~ from start and end
  })
  
  dist <- adist(maskedGeneIds$drerio, maskedGeneIds)
  if(output=='distance'){
    return(dist)
  }
  sd <- sd(dist)
  if(output=='sd'){
    return(sd)
  }
  gps <- 100/((sd+1)*fac)
  return(gps)
}

