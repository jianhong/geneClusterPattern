#' Drop non-standard chromosomes
#' @param genesList A list of GRanges object
#' @param negativePattern The pattern to match the seqnames that does not need.
#' @return A list of filtered GRanges
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels
#' @export
#' @examples
#' library(GenomicRanges)
#' gr <- GRanges(rep(c("chr2", "chr3", "chrM"), 2), IRanges(1:6, 10))
#' pruningSequences(list(gr))
#'
pruningSequences <- function(genesList, negativePattern='[_.M]'){
  if(is.list(genesList)){
    lapply(genesList, function(.ele){
      keepSeqlevels(.ele,
                    value=seqlevels(.ele)[
                      !grepl(negativePattern, seqlevels(.ele))],
                    pruning.mode = 'coarse')
    })
  }else{
    keepSeqlevels(genesList,
                  value=seqlevels(genesList)[
                    !grepl(negativePattern, seqlevels(genesList))],
                  pruning.mode = 'coarse')
  }
}
