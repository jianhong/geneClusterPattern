#' Apply colors to gene
#' @noRd
#' @param gr A GRanges object
#' @param colors The color sets
#' @return A GRanges object with colors
#' @importFrom BiocGenerics start end width strand
#' @importFrom grDevices palette
#' @importFrom stats setNames
addGeneInfo <- function(gr, colors){
  stopifnot(is(gr, 'GRanges'))
  stopifnot(length(names(gr))==length(gr))
  if(is.numeric(colors)){
    stopifnot(all(colors<9))
    ## convert number color to character
    colors <- setNames(palette()[colors], names(colors))
  }
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
