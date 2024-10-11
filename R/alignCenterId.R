#' @noRd
#' @param geneClusterPatterns gene patterns list
#' @param patternRegions rescaled patternRegions
#' @param id0 center id
#' @importFrom IRanges shift
#' @importFrom BiocGenerics width<-
#' @return A list, first element should be patterns, second element is
#' regions
alignCenterId <- function(geneClusterPatterns, patternRegions, id0){
  if(missing(id0)){
    genes <- lapply(geneClusterPatterns, names)
    gene_table <- table(unlist(genes))
    anchors <- names(gene_table)[gene_table==max(gene_table)]
    ## get percentage of center position of the anchors
    p0 <- vapply(geneClusterPatterns, function(pattern){
      if(any(anchors %in% names(pattern))){
        pattern <- pattern[names(pattern) %in% anchors]
        (min(start(pattern)) + max(end(pattern)))/2
      }else{
        0
      }
    }, numeric(1L))
  }else{
    ## get percentage of center position 
    p0 <- vapply(geneClusterPatterns, function(pattern){
      if(id0 %in% names(pattern)){
        (start(pattern[id0]) + end(pattern[id0]))/2
      }else{
        0
      }
    }, numeric(1L))
  }
  wid <- width(patternRegions)
  ## make the start position ratio same
  ## add the offset at the beginning and end
  ## (p0 + x)/(wid + x + y) = 0.5
  ## 2 * (p0 + x) = wid + x + y
  ## x - y = wid - 2*p0
  ## if (wid - 2*p0) < 0, set x=0, y=2*p0-wid
  ## if (wid - 2*p0) >= 0, set y=0, x=wid - 2*p0
  m <- wid - 2*p0
  x <- ifelse(m<0, 0, floor(m))
  y <- ifelse(m<0, ceiling(abs(m)), 0)
  geneClusterPatterns <- mapply(function(pattern, offset){
    unique(shift(pattern, shift = offset))
  }, geneClusterPatterns, x, SIMPLIFY = FALSE)
  width(patternRegions) <- width(patternRegions) + x + y
  return(list(
    patterns=geneClusterPatterns,
    regions=patternRegions
  ))
}
