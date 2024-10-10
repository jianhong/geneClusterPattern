#' retreive a homologs by given ensembl gene ids
#' @description
#' This function retrieves the user specified homologs for given ensembl gene
#' ids.
#' @param species The target species for homologs query.
#' @param mart object of class Mart.
#' @param ensembl_gene_ids Ensembl gene ids.
#' @param ... Parameters could be used by \link[biomaRt]{useEnsembl} except
#' `biomart`.
#' @return An object of GRanges
#' @importFrom biomaRt getBM useEnsembl
#' @export
#' @examples
#' if(interactive()){
#'   library(biomaRt)
#'   mart <- useMart('ensembl', 'drerio_gene_ensembl')
#'   homologsFromEnsemblIDs(
#'     'human',
#'     mart,
#'    c('ENSDARG00000008803', 'ENSDARG00000017038'))
#' }
#'
homologsFromEnsemblIDs <- function(species, mart, ensembl_gene_ids, ...){
  stopifnot(is.character(ensembl_gene_ids))
  stopifnot(is(mart, 'Mart'))
  species <- guessSpecies(species, output='abbr', ...)
  attributes <- paste0(species, '_homolog_ensembl_gene')
  bm <- getBM(attributes = c('ensembl_gene_id', attributes),
              filters = 'ensembl_gene_id',
              values = ensembl_gene_ids,
              mart = mart)
  bm <- bm[!is.na(bm[, 2]) & bm[, 2]!='', , drop=FALSE]
  mart2 <- useEnsembl('ensembl', paste0(species, '_gene_ensembl'), ...)
  gr <- grangesFromEnsemblIDs(mart2, ensembl_gene_ids=bm[, 2])
  if(length(gr)){
    gr$homolog_ensembl_gene_ids <- bm[match(names(gr), bm[, 2]),
                                      'ensembl_gene_id']
  }
  return(gr)
}
