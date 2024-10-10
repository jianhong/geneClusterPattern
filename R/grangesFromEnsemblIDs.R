#' create a GRanges object by given ensembl gene ids
#' @description
#' Create a GRanges object by given ensembl gene ids
#' @param mart object of class Mart.
#' @param ensembl_gene_ids Ensembl gene ids.
#' @return An object of GRanges
#' @importFrom biomaRt getBM
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
#' @examples
#' if(interactive()){
#'   library(biomaRt)
#'   mart <- useMart('ensembl', 'drerio_gene_ensembl')
#'   grangesFromEnsemblIDs(mart, c('ENSDARG00000008803', 'ENSDARG00000017038'))
#' }
#'
grangesFromEnsemblIDs <- function(mart, ensembl_gene_ids){
  stopifnot(is.character(ensembl_gene_ids))
  stopifnot(is(mart, 'Mart'))
  attributes <- c('ensembl_gene_id',
                  'external_gene_name',
                  'chromosome_name',
                  'start_position',
                  'end_position',
                  'strand')
  bm <- getBM(attributes = attributes,
              filters = 'ensembl_gene_id',
              values = ensembl_gene_ids,
              mart = mart)
  if(nrow(bm)>0){
    id <- as.character(bm[, 'ensembl_gene_id', drop=TRUE])
    chr <- as.character(bm[, 'chromosome_name', drop=TRUE])
    start <- as.numeric(bm[, 'start_position', drop=TRUE])
    end <- as.numeric(bm[, 'end_position', drop=TRUE])
    strand <- as.character(bm[, 'strand', drop=TRUE])
    strand <- ifelse(strand=='-1', '-', '+')
    keep <- !is.na(start)
    gr <- GRanges(seqnames = chr[keep],
                  ranges = IRanges(start = start[keep],
                                   end = end[keep],
                                   names = id[keep]),
                  strand = strand[keep],
                  gene_name = bm[keep, 'external_gene_name', drop=TRUE])
  }else{
    gr <- GRanges()
  }
  return(gr)
}
