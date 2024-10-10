#' @noRd
#' @param gr A GRanges object
#' @param colors The color sets
#' @return A GRanges object with colors
#' @importFrom BiocGenerics start end width strand
addGeneInfo <- function(gr, colors){
  gr$feature <- rep('exon', length(gr))
  gr$featureID <- gr$gene_name
  gr$hide_label <- !names(gr) %in% names(colors)
  gr$must_have_label <- !gr$hide_label
  gr$color <- ifelse(gr$must_have_label,
                     colors[names(gr)], '#33333311')
  grs <- lengths(lapply(split(strand(gr), gr$featureID), unique))
  if(length(grs)>0){
    ## remove the genes in multiple strand
    gr <- gr[gr$featureID %in% names(grs)[grs==1]]
  }
  return(gr)
}
