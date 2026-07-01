#' calculate the gene cluster
#' @description
#' A gene cluster is defined as the k nearest genes surrounding a reference gene
#' according to the rank matrix. The rank matrix is a gene-by-species matrix,
#' where each column contains the genomic ranks of genes for a given species.
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
#' fish <- readRDS(system.file('extdata', 'fish.rds',
#'                             package = 'geneClusterPattern'))
#' homologs <- readRDS(system.file('extdata', 'homologs.rds',
#'                                 package = 'geneClusterPattern'))
#' queryGene <- 'inhbaa'
#' nearestNeighbors <- getGeneCluster(fish, queryGene, homologs)
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
  # order the query genes
  queryRnk <- getGeneRank(queryGR)
  # order the homologs
  homoRnk <- lapply(homologsList, getGeneRank)
  # extract the shared chromosomes
  chr <- vapply(queryRnk, function(.ele){
    queryGeneName %in% names(.ele)
  }, logical(1L))
  if(any(chr)){
    if(sum(chr)==1){# keep the chromosome with the queryGeneName
      queryRnk <- queryRnk[[which(chr)]]
    }else{
      stop(queryGeneName, ' is listed in multiple chromosomes!')
    }
  }else{
    stop(queryGeneName, ' is not in names of queryRnk. Please report this bug!')
  }
  # get the rank position of the query gene
  center <- queryRnk[queryGeneName]
  # get the relative rank positions of query
  queryRnk0 <- queryRnk - center
  y0 <- lapply(homoRnk, function(y, center){
    # which chromosome contain the homolog
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
  # get all symbols involved, including query and homologs
  syms <- sort(unique(c(unlist(lapply(y0, names)), names(queryRnk))))
  # find the gene rank position in query
  queryRnk0 <- queryRnk0[syms]
  # find the gene rank positions in subjects
  y0 <- lapply(y0, function(.ele) .ele[syms])
  # create a matrix with all rank positions
  y <- do.call(cbind, c(list(queryRnk0=queryRnk0), y0))
  rownames(y) <- syms
  # make all un-available homologs as far distance
  y[is.na(y)] <- -length(syms)

  if(radius>0){
    k <- nrow(y)
    nearest <- nn2(y, k=k, radius = radius, searchtype = 'radius')
  }else{
    k <- min(k, nrow(y))
    nearest <- nn2(y, k=k)
  }
  # extract the nearest gene names
  ids <- rownames(y)[nearest$nn.idx[which(rownames(y)==names(center)), ,
                                    drop=TRUE]]
  ## make sure the query Id is the first one
  ids <- unique(c(names(center), ids))
  return(ids)
}
