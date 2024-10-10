#' @noRd
#' @param gr A GRanges
#' @return A GRanges with width rescaled.
#' @importFrom BiocGenerics start end width strand
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges ranges<-
rescalRegion <- function(gr){
  if(is.list(gr)) gr <- unlist(GRangesList(gr))
  width <- width(gr)
  width[width<1] <- 1
  scale <- max(width)/width
  width <- ceiling(width(gr)*scale/2)
  center <- ceiling((start(gr) + end(gr))/2)
  ranges(gr) <- IRanges(start = ifelse(center>width, center-width, 1),
                        end = center+width,
                        names = names(gr))
  return(gr)
}
