#' @noRd
#' @param gr A GRanges object.
#' @param ids IDs to be plotted. See \link{getGeneCluster}.
#' @param additionalID Other IDs should be included.
#' @param max_gap The maximal gaps from the first ID in `ids`.
#' @return a list with geneModel GRanges, region for that species,
#'  id0_width_pct and id_width_pct.
#' @importFrom BiocGenerics start<- end<-
#' @importFrom IRanges subsetByOverlaps
getGeneModel <- function(gr, ids, max_gap){
  stopifnot(length(ids)>1)
  id0 <- ids[1]
  w1 <- gr[ids[ids %in% names(gr)]]
  if(length(w1)==0){
    return(list(geneModel=GRanges('NA', IRanges(1, 1)),
                region=GRanges('NA', IRanges(1, 5000000)),
                id0_width_pct = 0,
                id_width_pct = 0)
    )
  }
  ## subset the gr by ids and max_gap
  strand(w1) <- '*'
  inRg <- id0 %in% names(gr)
  seqn <- ifelse(inRg, seqnames(gr[id0]), seqnames(w1)[1])
  w1 <- w1[as.character(seqnames(w1)) %in% as.character(seqn)]
  w2 <- range(w1)
  if(inRg) {
    w <- c(start(w2), end(w2)) - start(gr[id0])
    w <- max(abs(w))
    start(w2) <- start(gr[id0]) - w
    end(w2) <- end(gr[id0]) + w
  }
  region <- range(w2)
  if(inRg && !missing(max_gap)){
    cen <- gr[id0]
    start(region) <- max(1, start(cen)-max_gap, start(region))
    end(region) <- min(end(cen)+max_gap, end(region))
  }
  ext <-  min(round(width(region)/10), 5000)
  start(region) <- max(1, start(region) - ext)
  end(region) <- end(region) + ext
  w2 <- subsetByOverlaps(gr, w2)

  return(list(geneModel=w2, region=region,
              id0_width_pct = ifelse(
                id0 %in% names(w2),
                width(w2[id0])/width(region),
                0),
              id_width_pct = ifelse(
                any(ids %in% names(w2)),
                width(w2[names(w2) %in% ids])/width(region),
                0)
  )
  )
}
