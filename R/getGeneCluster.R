#' calculate the gene cluster
#' @description
#' A short description...
#' @param queryGR The query GRanges outputted by \link{grangesFromEnsemblIDs}.
#' @param queryGeneName The center gene used to search.
#' @param homologsList The homologs list outputted by \link{getHomologGeneList}.
#' @param k The maximum number of nearest genes.
#' @param radius The maximum number of gap genes for the search. If
#' radius is set to greater than 0, k parameter will be ignored.
#' @return The ids of the cluster for the center gene.
#' @importFrom stats dist
#' @importFrom RANN nn2
#' @importFrom methods is
#' @export
#' @examples
#' # example code
#'
getGeneCluster <- function(queryGR, queryGeneName,
                           homologsList, k = 10, radius = 0){
  stopifnot(length(queryGeneName)==1)
  stopifnot(is.list(homologsList))
  stopifnot(is.numeric(k))
  stopifnot(is.numeric(radius))
  stopifnot(length(k)==1)
  stopifnot(length(radius)==1)
  stopifnot(is(queryGR, 'GRanges'))
  if(!queryGeneName %in% names(queryGR)){
    stopifnot('gene_name should be a column of metadata for queryGR'=
                length(queryGR$gene_name)==length(queryGR))
    stopifnot('Could not find queryGeneName in queryGR'=
                queryGeneName %in% queryGR$gene_name)
    queryGeneName <- names(queryGR)[which(queryGR$gene_name==queryGeneName)]
  }
  queryRnk <- getGeneRank(queryGR)
  homoRnk <- lapply(homologsList, getGeneRank)
  chr <- vapply(queryRnk, function(.ele){
    queryGeneName %in% names(.ele)
  }, logical(1L))
  if(any(chr)){
    if(sum(chr)==1){
      queryRnk <- queryRnk[[which(chr)]]
    }else{
      stop(queryGeneName, ' is listed in multiple chromosomes!')
    }
  }else{
    stop(queryGeneName, ' is not in names of queryRnk. Please report this bug!')
  }

  center <- queryRnk[queryGeneName]
  queryRnk0 <- queryRnk - center
  y0 <- lapply(homoRnk, function(y, center){
    y <- y[vapply(y, function(.ele) names(center) %in% names(.ele),
                  FUN.VALUE = logical(1L))]
    if(length(y)>=1){
      ## find the chromosome with maximal genes and query gene
      y <- y[lengths(y)==max(lengths(y))][[1]]
      ## recenter the gene rank
      y0 <- y - y[names(center)]
    }else{
      y0 <- NA
    }
    return(y0)
  }, center=center)
  syms <- sort(unique(c(unlist(lapply(y0, names)), names(queryRnk))))
  queryRnk0 <- queryRnk0[syms]
  y0 <- lapply(y0, function(.ele) .ele[syms])
  y <- do.call(cbind, c(list(queryRnk0=queryRnk0), y0))
  rownames(y) <- syms
  y[is.na(y)] <- -length(syms)
  d <- dist(y)
  if(radius>0){
    k <- nrow(y)
    nearest <- nn2(y, k=k, radius = radius, searchtype = 'radius')
  }else{
    k <- min(k, nrow(y))
    nearest <- nn2(y, k=k)
  }
  ids <- rownames(y)[nearest$nn.idx[which(rownames(y)==names(center)), ,
                                    drop=TRUE]]
  ## make sure the query Id is the first one
  ids <- unique(c(names(center), ids))
  return(ids)
}
