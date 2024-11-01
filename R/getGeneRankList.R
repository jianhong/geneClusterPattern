#' Get the gene position orders for homologs
#' @description
#' Get the gene position orders without considering the strand info.
#' @param species The target species for homologs query.
#' @param mart object of class Mart.
#' @param ensembl_gene_ids Ensembl gene ids.
#' @param ... Parameters could be used by \link[biomaRt]{useEnsembl} except
#' `biomart`.
#' @return A list of list with gene ranks named by the chromosome names.
#' @export
#' @examples
#' # example code
#'
getHomologGeneRankList <- function(species, mart, ensembl_gene_ids, ...){
  homologs <- lapply(species, homologsFromEnsemblIDs, mart=mart,
                     ensembl_gene_ids = ensembl_gene_ids, ...)
  names(homologs) <- species
  ## convert the name and homolog_ensembl_gene_ids
  homologs <- convertNamesOfHomologIDs(homologs)
  
  homoRnk <- lapply(homologs, getGeneRank)
  return(homoRnk)
}

convertNamesOfHomologIDs <- function(homologs){
  lapply(homologs, function(.ele){
    n <- names(.ele)
    names(.ele) <- .ele$homolog_ensembl_gene_ids
    .ele$homolog_ensembl_gene_ids <- n
    .ele
  })
}