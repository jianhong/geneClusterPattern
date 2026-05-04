#' Read Orthologues results as a list
#' This function is under-development.
#' @param orthogroups A dataframe with the gene ids for each orthologues group
#' @param annoGR_list A GRangesList object with the gene information. 
#' The names of each element should be contain all the gene ids saved in the 
#' orthogroups parameter. For each element, the 'gene_name' should be available
#' in the metadata.
#' @param query_species The species used as queryGR in \link{getGeneCluster}.
#' @param split character vector used to split the gene ids in the dataframe of 
#' orthogroups parameter.
#' @return a list of queryGR and homologsList used by \link{getGeneCluster}.
#' @export
#' @examples
#' orthogroups <- read.delim(
#'      system.file('extdata/orthofinder/Orthogroups.tsv.gz',
#'                   package='geneClusterPattern'),
#'      row.names=1)
#' annoGR_list <- readRDS(
#'      system.file('extdata/orthofinder/grange.obj.rds',
#'                   package='geneClusterPattern'))
#' colnames(orthogroups) <- names(annoGR_list)
#' hres <- fromOrthogroups(orthogroups, annoGR_list)
#' 
fromOrthogroups <- 
  function(orthogroups, annoGR_list, query_species=1, split=', '){
    if(is.character(query_species)){
      stopifnot('query_species must in the names of annoGR_list'=
                  query_species %in% names(annoGR_list))
    }else{
      query_species <- names(annoGR_list)[query_species]
    }
    stopifnot(inherits(annoGR_list, c('GRangesList', 'list')))
    null <- lapply(annoGR_list, function(.ele){
      stopifnot('annoGR_list must be an object of GRangesList'=
                  is(.ele, 'GRanges'))
      stopifnot('annoGR_list elements must contain metadata "gene_name"'=
                  length(.ele$gene_name)==length(.ele))
      stopifnot('annoGR_list elements must contain names.'=
                  length(names(.ele))==length(.ele))
    })
    stopifnot(is.data.frame(orthogroups))
    stopifnot(identical(colnames(orthogroups), names(annoGR_list)))

  orthogroups.list <- apply(orthogroups, 2, function(.ele){
    strsplit(.ele, split)
  }, simplify = FALSE)
  
  sp1 <- unique(unlist(orthogroups.list[[query_species]], use.names = FALSE))
  trimName <- FALSE
  if(sum(sp1 %in% names(annoGR_list[[query_species]]))<0.5*length(sp1)){
    sp1 <- trimENSname(sp1)
    trimName <- TRUE
  }
  if(sum(sp1 %in% names(annoGR_list[[query_species]]))<0.1*length(sp1)){
    stop('less than 10% gene info for the query species exist in annoGR_list')
  }
  sp1 <- annoGR_list[[query_species]][sp1]
  
  query <- orthogroups.list[[query_species]]
  l <- lengths(query)
  homologs <- mapply(FUN=function(target, b){
    l1 <- lengths(target)
    q <- mapply(function(.e, .l) {
      rep(.e, .l)
    }, query, l1)
    p <- mapply(function(.e, .l){
      rep(.e, each=.l)
    }, target, l)
    stopifnot(identical(lengths(p), lengths(q)))
    q <- unlist(q, use.names = FALSE)
    p <- unlist(p, use.names = FALSE)
    if(trimName){
      q <- trimENSname(q)
      p <- trimENSname(p)
    }
    b <- b[p]
    b$homolog_ensembl_gene_ids <- q
    b
  }, orthogroups.list[!names(annoGR_list) %in% query_species],
  annoGR_list[!names(annoGR_list) %in% query_species], SIMPLIFY = FALSE)
  homologs <- convertNamesOfHomologIDs(homologs)
  return(
    list(queryGR=sp1, homologsList=homologs)
  )
}

trimENSname <- function(x) sub('\\.\\d+$', '', x)





