#' Calculate the gene order score
#' @description
#' The gene order score is the mean of 
#' edit (Levenshtein) distance of given ids from different species.
#' Gene order score = mean(adist)*k/length(ids).
#' The higher the gene order score, the lower the conservation of the gene order
#' @param genesList A list of GRanges.
#' @param ids IDs to be plotted. See \link{getGeneCluster}.
#' @param ref The reference species.
#' @param k The maximum number of nearest genes used in \link{getGeneCluster}.
#' @param max_gap The maximal gaps from the first ID in `ids`.
#' @param output Output values
#' @return the standard deviation or edit distance values.
#' @importFrom stats sd
#' @export
#' @examples
#' # example code
#' 
geneOrderScore <- function(genesList, ids, ref, k=length(ids), max_gap=1e7,
                             output=c("score", "distance")){
  output <- match.arg(output)
  stopifnot(is.numeric(k))
  stopifnot(length(k)==1)
  if(!missing(ref)){
    stopifnot(ref %in% names(genesList))
  }
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
  
  # get the strings of gene order for each species
  geneModels <- checkAndGetGeneModels(genesList, ids=ids, max_gap=max_gap)
  geneIds <- lapply(geneModels, function(.ele) names(.ele$geneModel))
  
  # id map
  geneIdMap <- geneIdMap[seq_along(ids)]
  names(geneIdMap) <- ids
  
  ## mask the non-essential IDs
  stringList <- lapply(geneIds, function(.ele) geneIdMap[.ele])
  
  ## rever string if reversed
  stringList <- reverseByCor(stringList)
  
  maskedGeneIds <- lapply(stringList, function(ch){
    ch[is.na(ch)] <- '~' # 126
    sub('^~+', '',
        sub('~+$', '',
            paste(ch, collapse=''))) # remove ~ from start and end
  })
  
  
  if(missing(ref)){
    dist <- adist(maskedGeneIds)
  }else{
    dist <- adist(maskedGeneIds[[ref]], maskedGeneIds)
    dist <- dist[, colnames(dist)!=ref, drop=TRUE]
  }
  if(output=='distance'){
    return(dist)
  }
  if(missing(ref)){
    gps <- mean(dist[lower.tri(dist, diag = FALSE)], na.rm=TRUE)*fac
  }else{
    gps <- mean(dist, na.rm=TRUE)*fac
  }
  return(gps)
}

reverseByCor <- function(stringList){
  a <- unique(unlist(stringList))
  a <- a[!is.na(a)]
  b <- lapply(stringList, function(.ele){
    match(a, .ele)
  })
  cor <- lapply(seq_along(b)[-1], function(.ele)
    tryCatch(cor(b[[1]], b[[.ele]], method = 'spearman'),
             error=function(.e){
               0
             }))
  if(sum(cor<0)>=length(cor)/2){
    ## reverse the strings with negative cor
    k <- which(cor<0)+1
  }else{
    ## reverse the strings with positive cor
    k <- c(1, which(cor>0)+1)
  }
  stringList[k] <- lapply(stringList[k], rev)
  return(stringList)
}
