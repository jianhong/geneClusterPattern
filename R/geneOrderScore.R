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
#' 'spearman correlation', the score is the mean of absolute value of
#' Spearman correlation of the genes appearance order.
#' 'non-random score', the score is the mean of Jaccard index of 2-order and
#' 3-order of ids for paired samples.
#' The higher the gene order score, the higher the conservation of the gene order.
#' @param output Output values
#' @return the score or values matrix such as adist, or alignment scores.
#' @importFrom stats sd
#' @importFrom pwalign pairwiseAlignment
#' @importFrom Biostrings AA_PROTEINOGENIC
#' @export
#' @examples
#' # example code
#' fish <- readRDS(system.file('extdata', 'fish.rds',
#' package = 'geneClusterPattern'))
#' homologs <- readRDS(system.file('extdata', 'homologs.rds',
#'                                 package = 'geneClusterPattern'))
#' queryGene <- 'pcolce2b'
#' nearest10neighbors <- getGeneCluster(fish, queryGene, homologs, k=10)
#' genesList <- c(drerio=fish, homologs)[
#' c("hsapiens", "mmusculus", "drerio", "olatipes", "nfurzeri", "gaculeatus")]
#' geneOrderScore(genesList, nearest10neighbors, ref='drerio')
#' geneOrderScore(genesList, nearest10neighbors, method='non-random score')
geneOrderScore <- function(genesList, ids, ref, k=length(ids), max_gap=1e7,
                           method = c('edit distance',
                                      'global alignment score',
                                      'spearman correlation',
                                      'non-random score'),
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
  # id map
  geneIdMap <- geneIdMap[seq_along(ids)]
  names(geneIdMap) <- ids
  
  # get the strings of gene order for each species
  geneModels <- checkAndGetGeneModels(genesList, ids=ids, max_gap=max_gap)
  geneIds <- lapply(geneModels, function(.ele) names(.ele$geneModel))
  #strand <- lapply(geneModels, function(.ele) strand(.ele$geneModel))

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
            paste(ch, collapse=''))) # remove gapLetter from start and end
  })
  
  FUN <- get(gsub(' |-', '_', method))
  dist <- FUN(stringList, maskedGeneIds, ref)
  
  if(output!='score'){
    return(dist)
  }
  
  if(method=='spearman correlation'){## convert to positive number
    dist <- abs(dist)
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
    dist <- unlist(dist)
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
    cor <- unlist(cor)
  }
  cor[is.na(cor)] <- 0
  return(cor)
}

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

non_random_score <- function(stringList, maskedGeneIds, ref){
  ## gene random distributed, 
  ## 2 order Markov chain, or 3 order Markov chain
  # get all pairs in a and b
  # check the shared pairs
  # index = shared pairs / all pairs
  # like Jaccard index
  b <- lapply(maskedGeneIds, function(.ele){
    .ele <- gsub('--+', '--', .ele) # gapLetter <- '-' # 45
    n <- seq.int(nchar(.ele))
    s <- c(strsplit(.ele, '')[[1]],
           substring(.ele, first =n[-length(n)], last=n[-1]))
    if(length(n)>2){
      s <- c(s, substring(.ele,
                          first =n[-c(length(n)-1, length(n))],
                          last=n[-c(1, 2)]))
    }
    s <- strsplit(s, '')
    s <- lapply(s, function(.e){
      paste(sort(unique(.e[.e!='-'])), collapse = '')
    })
    s <- unlist(unique(s))
    s[s!='']
  })
  if(missing(ref)){
    comb <- getComb(names(b))
    ji <- apply(comb, 1, function(.ele) 
      jaccard(b[[.ele[1]]], b[[.ele[2]]]))
    ji <- matrix(ji, nrow = length(b),
                  dimnames = list(names(b), names(b)))
  }else{
    ji <- lapply(seq_along(b), function(.ele)
      jaccard(b[[ref]], b[[.ele]]))
    ji <- unlist(ji)
  }
  return(ji)
}

kolmogorov_complexity <- function(stringList, maskedGeneIds, ref){
  
}