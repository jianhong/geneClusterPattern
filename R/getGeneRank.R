#' Get the gene position orders
#' @description
#' Get the gene position orders without considering the strand info.
#' @param genes A GRanges object with names.
#' @return A list with gene ranks named by the chromosome names.
#' @importFrom IRanges promoters
#' @importFrom BiocGenerics strand
#' @importFrom GenomeInfoDb seqnames
#' @export
#' @examples
#' # example code
#'
getGeneRank <- function(genes){
  stopifnot(is(genes, 'GRanges'))
  stopifnot('Input must have name.'=length(names(genes))==length(genes))
  strand(genes) <- '*'
  genes <- promoters(genes, upstream = 0, downstream = 1)
  genes <- split(genes, as.character(seqnames(genes)))
  genes <- lapply(genes, sort)
  genesRnk <- lapply(genes, function(.ele){
    rnk <- seq_along(.ele)
    names(rnk) <- names(.ele)
    return(rnk)
  })
  return(genesRnk)
}
