pasteReplaceLast <- function(..., sep=' ', collapse=', ', last=', or'){
  out <- paste(..., sep=sep, collapse=collapse)
  out <- sub(', ([^,]*)$', paste0(last, ' \\1'), out)
  return(out)
}