#' @noRd
#' @param gr A GRanges object
#' @param region The region to plot
#' @param species The y label
#' @param scaleBar Add scale bar or not
#' @return NULL
#' @importFrom trackViewer trackList optimizeStyle setTrackViewerStyleParam
#' viewTracks
#' @importFrom methods new
plotOneTrack <- function(gr, region, species, scaleBar=FALSE){
  ## remove the gr in multiple strand
  grs <- lengths(lapply(split(strand(gr), gr$featureID), unique))
  gr <- gr[gr$featureID %in% names(grs)[grs==1]]
  if(length(gr)==0){
    gr <- addGeneInfo(region, c())
  }
  if(length(gr)){
    gene <- new('track', dat=gr, type='gene', name=species[1],
                style=new('trackStyle', color='lightblue'))
    tr <- trackList(gene)
    names(tr) <- species[1]
    optSty <- optimizeStyle(tr)
    trackList <- optSty$tracks
    viewerStyle <- optSty$style
    if(scaleBar){
      setTrackViewerStyleParam(viewerStyle, "xaxis", TRUE)
      setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .01, .05))
    }else{
      setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
      setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .05, .01, .05))
    }
    viewTracks(trackList, viewerStyle = viewerStyle, gr=region, newpage = FALSE)
  }
  return(invisible(NULL))
}
