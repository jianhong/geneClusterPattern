#' Get the gene position orders for homologs
#' @description
#' Get the gene position orders without considering the strand info.
#' @param species The target species for homologs query.
#' @param mart object of class Mart.
#' @param ensembl_gene_ids Ensembl gene ids.
#' @param ... Parameters could be used by \link[biomaRt]{useEnsembl} except
#' `biomart`.
#' @return A list of list with gene ranks named by the chromosome names.
#' @noRd
#' @examples
#' if(interactive()){
#'   ## Ensembl server may not response
#'   library(biomaRt)
#'   fish <- readRDS(system.file('extdata', 'fish.rds',
#'                             package = 'geneClusterPattern'))
#'   ensembl_gene_ids <- names(fish[seqnames(fish)=='24'])
#'   species <- c('hsapiens', 'mmusculus')
#'   fish_mart <- useMart("ENSEMBL_MART_ENSEMBL", "drerio_gene_ensembl")
#'   homologs <- getHomologGeneRankList(species, fish_mart, ensembl_gene_ids)
#' }
#' 
getHomologGeneRankList <- function(species, mart, ensembl_gene_ids, ...){
  homologs <- getHomologGeneList(species, mart, ensembl_gene_ids, ...)
  
  homoRnk <- lapply(homologs, getGeneRank)
  return(homoRnk)
}

