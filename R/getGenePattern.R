#' @noRd
#' @param gr A GRanges object with colors, output from addGeneInfo
#' @param reg The regions to plot
#' @param k The maximal non-essential genes to be plotted between two
#'  essential genes.
#' @return gene patter in GRanges object
#' @importFrom IRanges gaps promoters ranges<- disjoin reduce ranges
#' @importFrom GenomicRanges GRangesList
#' @importFrom BiocGenerics start end width strand strand<-
getGeneClusterPattern <- function(gr, reg, k=5){
  if(length(gr)){
    ## use simple plot
    gr$feature[gr$must_have_label] <- 'gene'
    ## sample the gr if the not essential genes is more than 5
    gr <- gr[order(start(promoters(gr, upstream=0, downstream=1)))]
    hide_label <- rle(gr$hide_label)
    keep <- mapply(function(v, l){
      true <- rep(TRUE, l)
      if(v){
        if(l>k){
          false <- sample(seq.int(l), size=l-k)
          true[false] <- FALSE
        }
      }
      return(true)
    }, hide_label$values, hide_label$lengths, SIMPLIFY = FALSE)
    keep <- unlist(keep)
    gr <- gr[keep]
    ## rescale the genes and the gaps
    grnostrand <- gr
    strand(grnostrand) <- '*'
    grgap <- gaps(grnostrand, start = start(reg), end = end(reg))
    grgap <- grgap[as.character(strand(grgap))=='*' &
                     as.character(seqnames(grgap))==
                     as.character(seqnames(reg)[1])]
    grgap <- unique(grgap)
    names(grgap) <- paste0('gap_', seq_along(grgap))
    grgap$hide_label <- TRUE
    grgap$must_have_label <- !grgap$hide_label
    grs <- c(gr, grgap)
    grs <- grs[order(start(promoters(grs, upstream=0, downstream=1)))]
    grs_dj <- disjoin(grs, with.revmap=TRUE, ignore.strand=TRUE)
    grs_dj$must_have_label <- vapply(grs_dj$revmap, function(.ele){
      any(grs[.ele]$must_have_label)
    }, FUN.VALUE = logical(1L))
    wid <- ifelse(!grs_dj$must_have_label,
                  ceiling(log10(width(grs_dj))),# log10 re-scale for non-essential
                  ceiling(log2(width(grs_dj)))) # log2 re-scale for essential
    pos <- cumsum(c(4, wid))
    ranges(grs_dj) <- IRanges(start = pos[seq_along(grs_dj)],
                              end = pos[-1])
    ## split the disjoined granges by the reverse map
    grs_dj <- split(rep(grs_dj, lengths(grs_dj$revmap)), unlist(grs_dj$revmap))
    grs_dj <- GRangesList(grs_dj)
    grs_dj <- reduce(grs_dj)
    grs_dj <- unlist(grs_dj)
    names(grs_dj) <- names(grs)
    ranges(grs) <- ranges(grs_dj)
    ## get strand information back
    keep <- names(gr) %in% names(grs)
    gr <- gr[keep]
    strand(grs[names(gr)]) <- strand(gr)
    gr <- grs[names(gr)]
  }
  return(gr)
}
