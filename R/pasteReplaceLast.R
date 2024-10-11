pasteReplaceLast <- function(..., sep=' ', collapse=', ', last=', or'){
  out <- paste(..., sep=sep, collapse=collapse)
  out <- sub(', ([^,]*)$', ', or \\1', out)
  return(out)
}