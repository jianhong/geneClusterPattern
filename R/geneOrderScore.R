#' Calculate the gene order score
#' @description
#' Calculate the gene order score.
#' @param genesList A list of GRanges.
#' @param ids IDs to be plotted. See \link{getGeneCluster}.
#' @param ref The reference species.
#' @param k The maximum number of nearest genes used in \link{getGeneCluster}.
#' @param max_gap The maximal gaps from the first ID in `ids`.
#' @param method The method to calculate the gene order score.
#' 'edit distance', the gene order score is the mean of 
#' edit (Levenshtein) distance of given ids from different species.
#' Gene order score = mean(adist)*k/length(ids).
#' The higher the gene order score, the lower the conservation of the gene order.
#' 'global alignment score', the gene order score is the mean of global 
#' alignment score. see \link[pwalign]{pairwiseAlignment}.
#' The higher the gene order score, the higher the conservation of the gene order.
#' 'spearman correlation', the score is the mean of Spearman correlation of
#' the genes appearance order.
#' @param output Output values
#' @return the score or values matrix such as adist, or alignment scores.
#' @importFrom stats sd
#' @importFrom pwalign pairwiseAlignment
#' @importFrom Biostrings AA_PROTEINOGENIC
#' @export
#' @examples
#' # example code
#' 
geneOrderScore <- function(genesList, ids, ref, k=length(ids), max_gap=1e7,
                           method = c('edit distance',
                                      'global alignment score',
                                      'spearman correlation'),
                           output=c("score", "value matrix")){
  output <- match.arg(output)
  method <- match.arg(method)
  stopifnot(is.numeric(k))
  stopifnot(length(k)==1)
  if(!missing(ref)){
    stopifnot(ref %in% names(genesList))
  }
  fac <- k/length(ids)
  gapLetter <- '-' # 45
  if(method=='global alignment score'){
    geneIdMap <- AA_PROTEINOGENIC
  }else{
    geneIdMap <- c(seq(65, 90),# A-Z
                   seq(97, 122), # a-z
                   seq(48, 57), # 0-9
                   seq(32, 44),
                   c(46, 47),
                   seq(58, 64),
                   seq(91, 96),
                   seq(123, 126))
    geneIdMap <- intToUtf8(geneIdMap)
    geneIdMap <- strsplit(geneIdMap, '')[[1]]
  }
  
  
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
  stringList <-
    reverseByCor(stringList,
                 ref = if(missing(ref)) names(stringList)[1]
                       else ref)
  
  maskedGeneIds <- lapply(stringList, function(ch){
    ch[is.na(ch)] <- gapLetter # 45
    sub(paste0('^', gapLetter, '+'), '',
        sub(paste0(gapLetter, '+$'), '',
            paste(ch, collapse=''))) # remove ~ from start and end
  })
  
  FUN <- get(gsub(' ', '_', method))
  dist <- FUN(stringList, maskedGeneIds, ref)
  
  if(output!='score'){
    return(dist)
  }
  if(is.matrix(dist)){
    gps <- mean(dist[lower.tri(dist, diag = FALSE)], na.rm=TRUE)*fac
  }else{
    gps <- mean(dist, na.rm=TRUE)*fac
  }
  return(gps)
}

reverseByCor <- function(stringList, ref){
  cor <- spearman_correlation(stringList, NULL, ref)
  if(sum(cor<0)>=length(cor)/2){
    ## reverse the strings with negative cor
    k <- which(cor<0)
  }else{
    ## reverse the strings with positive cor
    k <- which(cor>0)
  }
  stringList[k] <- lapply(stringList[k], rev)
  return(stringList)
}

edit_distance <- function(stringList, maskedGeneIds, ref){
  if(missing(ref)){
    dist <- adist(maskedGeneIds, ignore.case = FALSE, partial = FALSE)
  }else{
    dist <- adist(maskedGeneIds[[ref]], maskedGeneIds,
                  ignore.case = FALSE, partial = FALSE)
    dist <- dist[, colnames(dist)!=ref, drop=TRUE]
  }
  return(dist)
}

getComb <- function(x){
  expand.grid(x, x)
}

global_alignment_score <- function(stringList, maskedGeneIds, ref){
  if(missing(ref)){
    comb <- getComb(names(maskedGeneIds))
    dist <- apply(comb, 1, function(.ele){
      pairwiseAlignment(maskedGeneIds[[.ele[1]]],
                        maskedGeneIds[[.ele[2]]],
                        scoreOnly = TRUE, type='global')
    })
    dist <- matrix(dist, nrow = length(maskedGeneIds),
                   dimnames = list(names(maskedGeneIds), names(maskedGeneIds)))
  }else{
    dist <- lapply(maskedGeneIds[names(maskedGeneIds)!=ref], function(s2){
      pairwiseAlignment(maskedGeneIds[[ref]], s2, scoreOnly = TRUE, type='global')
    })
  }
  return(dist)
}
spearman_correlation <- function(stringList, maskedGeneIds, ref){
  a <- unique(unlist(stringList))
  a <- a[!is.na(a)]
  b <- lapply(stringList, function(.ele){
    match(a, .ele)
  })
  b <- lapply(b, function(.ele){
    .ele[is.na(.ele)] <- 0
    .ele
  })
  if(missing(ref)){
    comb <- getComb(names(b))
    cor <- apply(comb, 1, function(.ele){
      tryCatch(cor(b[[.ele[1]]], b[[.ele[2]]], method = 'spearman'),
               error=function(.e){
                 0
               })
    })
    cor <- matrix(cor, nrow = length(b),
                  dimnames = list(names(b), names(b)))
  }else{
    cor <- lapply(seq_along(b), function(.ele)
      tryCatch(cor(b[[ref]], b[[.ele]], method = 'spearman'),
               error=function(.e){
                 0
               }))
  }
  cor <- unlist(cor)
  cor[is.na(cor)] <- 0
  return(cor)
}
