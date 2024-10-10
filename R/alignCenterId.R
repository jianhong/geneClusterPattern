#' @noRd
#' @param geneClusterPatterns gene patterns list
#' @param patternRegions rescaled patternRegions
#' @param id0 center id
#' @importFrom IRanges shift
#' @importFrom BiocGenerics width<-
#' @return A list, first element should be patterns, second element is
#' regions
alignCenterId <- function(geneClusterPatterns, patternRegions, id0){
  p0 <- vapply(geneClusterPatterns, function(pattern){
    if(id0 %in% names(pattern)){
      start(pattern[id0])
    }else{
      0
    }
  }, numeric(1L))
  m <- max(p0)
  p1 <- ifelse(p0==0, 0, m - p0)
  geneClusterPatterns <- mapply(function(pattern, offset){
    unique(shift(pattern, shift = offset))
  }, geneClusterPatterns, p1, SIMPLIFY = FALSE)
  width(patternRegions) <- width(patternRegions) + m
  return(list(
    patterns=geneClusterPatterns,
    regions=patternRegions
  ))
}
