#' plot gene patters for given ids
#' @description
#' A short description...
#' @param genesList A list of GRanges.
#' @param ids IDs to be plotted. See \link{getGeneCluster}.
#' @param additionalID Other IDs should be included.
#' @param max_gap The maximal gaps from the first ID in `ids`.
#' @param colors The colors for genes
#' @param maxNonEssential The maximal number for non-essential genes between
#' two essential genes.
#' @return invisible list of plot data.
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom grid grid.newpage viewport pushViewport popViewport grid.text
#' @importFrom grDevices rainbow
#' @importFrom IRanges subsetByOverlaps ranges<-
#' @export
#' @examples
#' # example code
#'

plotGeneClusterPatterns <- function(genesList, ids, additionalID, max_gap=1e7,
                             colors, maxNonEssential=5){
  if(!missing(additionalID)) ids <- unique(c(ids, additionalID))
  if(missing(colors)){
    colors <- rainbow(length(ids))
    names(colors) <- ids
  }else{
    stopifnot('Not all ids are in names of colors'=
                all(ids %in% names(colors)))
  }
  stopifnot('genesList must have names'=
              length(names(genesList))==length(genesList))
  null <- lapply(genesList, function(.ele){
    if(!is(.ele, 'GRanges')){
      stop('genesList must be a list of GRanges', call. = FALSE)
    }
    if(length(.ele$gene_name)!=length(.ele)){
      stop('Elements in genesList must have metadata "gene_name"',
           call. = FALSE)
    }
    if(length(names(.ele))!=length(.ele)){
      stop('Elements in genesList must have names',
           call. = FALSE)
    }
  })

  geneModels <- lapply(genesList, getGeneModel,
                       ids=ids, max_gap=max_gap)
  ## extract the plot regions for each species
  region <- unlist(GRangesList(lapply(geneModels, function(.ele) .ele$region)))
  ## make the plot region size comparable
  region <- rescalRegion(region)
  ## extract the genes in the rescaled regions
  geneModels <- mapply(function(.ele, .region){
    w2 <- subsetByOverlaps(.ele, .region)
    if(length(w2)>2*length(ids)){
      n2 <- names(w2)[!names(w2) %in% ids]
      n2 <- sample(n2, length(ids))
      w2 <- w2[names(w2) %in% c(ids, n2)]
    }
    w2
  }, genesList, split(region, names(region))[names(genesList)])
  ## apply colors
  geneModels <- lapply(geneModels, addGeneInfo, colors=colors)
  ## create the gene patterns
  geneClusterPatterns <- mapply(getGeneClusterPattern, geneModels,
                         split(region, names(region))[names(geneModels)],
                         k=maxNonEssential)
  ## get the plot region for gene patterns
  patternRegions <- lapply(geneClusterPatterns, function(gr){
    if(length(gr)){
      strand(gr) <- '*'
      reg <- range(gr)
      ranges(reg) <- IRanges(start = 1, end = end(reg)+3)
      return(reg)
    }else{
      return(GRanges('NA', IRanges(1, 50)))
    }
  })
  ## make the plot region size comparable
  patternRegions <- rescalRegion(patternRegions)
  ## align center
  geneClusterPatterns <- alignCenterId(geneClusterPatterns, patternRegions, ids[1])
  patternRegions <- geneClusterPatterns$regions
  geneClusterPatterns <- geneClusterPatterns$patterns
  h <- 1/length(genesList)
  grid.newpage()
  grid.text('Gene Patterns', x = 0.5, y=0.975)
  pushViewport(viewport(y=0.525, height = 0.45, just=c(0.5, 0)))
  for(i in seq_along(geneClusterPatterns)){
    vp <- viewport(x=0, y=h*(i-1), width=0.95, height = h, just = c(0, 0))
    pushViewport(vp)
    plotOneTrack(geneClusterPatterns[[i]],
              region=patternRegions[i],
              species=names(geneClusterPatterns)[i])
    popViewport()## vp
  }
  popViewport()
  grid.text('Real Gene Track', x = 0.5, y=0.5)
  pushViewport(viewport(y=0.05, height = 0.45, just=c(0.5, 0)))
  for(i in seq_along(geneModels)){
    vp <- viewport(x=0, y=h*(i-1), width=0.95, height = h, just = c(0, 0))
    pushViewport(vp)
    plotOneTrack(geneModels[[i]],
              region=region[names(geneModels)[i]],
              species=paste0(names(geneModels)[i], ': chr',
                             as.character(seqnames(region[names(geneModels)[i]])))[1],
              scaleBar = TRUE)
    popViewport()## vp
  }
  popViewport()
  return(invisible(list(
    geneClusterPatterns=geneClusterPatterns,
    patternRegions=patternRegions,
    geneModels=geneModels,
    geneRegions=region
  )))
}
