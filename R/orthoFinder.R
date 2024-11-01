#' Read orthoFinder Orthologues results as a list
#' @param path The path of the orthologues result file eg:
#' 'primary_transcripts/OrthoFinder/Results_Nov01/Orthologues/Danio_rerio.GRCz11.pep.all.tsv'
#' @return a list of matrix with the paired orthologues.
#' @importFrom utils read.delim
#' @export  
orthologPairsFromOrthoFinder <- function(path){
  stopifnot(file.exists(path))
  x <- read.delim(path)
  x <- split(x, x$Species)
  x <- lapply(x, function(.ele){
    cn <- c(colnames(.ele)[3], .ele[1, 2])
    a <- strsplit(.ele[, 3], ', ')
    b <- strsplit(.ele[, 4], ', ')
    .ele <- mapply(a, b, FUN=function(.a, .b){
      cbind(x=rep(.a, each=length(.b)), y=rep(.b, length(.a)))
    })
    .ele <- unique(do.call(rbind, .ele))
    colnames(.ele) <- cn
    .ele
  })
  return(x)
}

#' Get the gene positions for output of OrthoFinder
#' @description
#' Get the gene positions without considering the strand info.
#' @param orthologs a two columns matrix
#' @param mart An Mart object
#' @return A GRanges object with homolog_ensembl_gene_ids
#' @export
#' @examples
#' # example code
#' 
getHomologForOrthoFinder <- function(orthologs, mart){
  stopifnot(is(mart, 'Mart'))
  stopifnot(
    'orthologs must be an element of output of orthorlogPairsFromOrthoFinder'=
              is.matrix(orthologs))
  stopifnot(
    'orthologs must be an element of output of orthorlogPairsFromOrthoFinder'=
              ncol(orthologs)==2)
  orthologs <- sub('\\.\\d+$', '', orthologs)
  stopifnot('The gene id is not ensembl ids'=
              all(grepl('ENS', orthologs[, 2])))
  gr <- grangesFromEnsemblIDs(mart = mart, ensembl_gene_ids = orthologs[, 2])
  if(length(gr)){
    gr$homolog_ensembl_gene_ids <- 
      orthologs[match(names(gr), orthologs[, 2]), 1]
  }
  return(gr)
}


#' Get the gene position orders for output of OrthoFinder
#' @description
#' Get the gene position orders without considering the strand info.
#' @param orthologs The output of \link{orthologPairsFromOrthoFinder}.
#' @param marts a list of object of class Mart.
#' @return A list of list with gene ranks named by the chromosome names.
#' @export
#' @examples
#' # example code
#'
getHomologListForOrthoFinder <- function(orthologs, marts){
  stopifnot(
    'orthologs must be an element of output of orthorlogPairsFromOrthoFinder'=
      is.list(orthologs))
  stopifnot(length(marts)==length(orthologs))
  homologs <- mapply(getHomologForOrthoFinder, orthologs, marts)
  ## convert the name and homolog_ensembl_gene_ids
  homologs <- convertNamesOfHomologIDs(homologs)
  return(homologs)
}