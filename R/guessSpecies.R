#' Guess the species
#' @description
#' Guess the species for a given string.
#' @param species A character(1L) to guess.
#' @param output Output type.
#' @param ... Parameters could be used by \link[biomaRt]{useEnsembl} except
#' `biomart`.
#' @return A character or Mart object.
#' @importFrom methods getPackageName
#' @importFrom biomaRt useEnsembl listDatasets
#' @importFrom utils adist
#' @export
#' @examples
#' if(interactive()){
#'   guessSpecies('stickleback')
#' }
#'
guessSpecies <- function(species,
                         output=c('abbr', 'scientific name',
                                  'mart',
                                  'taxid', 'common name'),
                         ...){
  stopifnot(is.character(species))
  stopifnot(length(species)==1)
  output <- match.arg(output)
  mart <- useEnsembl(biomart = 'ensembl', ...)
  dataset <- listDatasets(mart = mart)
  taxonomy <- readRDS(system.file('extdata', 'taxonomy_names.rds',
                                  package = getPackageName()))
  abbr <- taxonomy$scientific_name
  abbr <- do.call(rbind, strsplit(abbr, ' '))
  taxonomy$abbr <- tolower(paste0(substr(abbr[, 1], 1, 1), abbr[, 2]))
  if(paste0(species, '_gene_ensembl') %in% dataset$dataset){
    id <- which(taxonomy$abbr==species & !is.na(taxonomy$abbr))
    return(switch(output,
                  'abbr'=species,
                  'mart'=useEnsembl(biomart = 'ensembl',
                                    dataset=paste0(species, '_gene_ensembl')),
                  'scientific name'=taxonomy[id[1], 'scientific_name'],
                  'taxid'=taxonomy[id[1], 'taxid'],
                  'common name'=taxonomy[id[1], 'common_name']))
  }
  dist_scientific_name <- adist(species, taxonomy$scientific_name,
                                ignore.case = TRUE)[1, ]
  dist_common_name <- adist(species, taxonomy$common_name,
                            ignore.case = TRUE)[1, ]
  min_scientific_name <- min(dist_scientific_name, na.rm = TRUE)
  min_common_name <- min(dist_common_name, na.rm = TRUE)
  if(min_scientific_name==0){
    id <- which(dist_scientific_name==min_scientific_name &
                  !is.na(dist_scientific_name))
  }else if(min_common_name==0){
    id <- which(dist_common_name==min_common_name & !is.na(dist_common_name))
  }else if(min_scientific_name<min_common_name){
    id <- which(dist_scientific_name==min_scientific_name &
                  !is.na(dist_scientific_name))
    warning('Can not find exactly match. Using ',
            taxonomy[id[1], 'scientific_name'])
  }else{
    id <- which(dist_common_name==min_common_name & !is.na(dist_common_name))
    warning('Can not find exactly match. Using ',
            taxonomy[id[1], 'common_name'])
  }
  sname0 <- taxonomy[id[1], 'scientific_name']
  sname <- strsplit(sname0, ' ')[[1]]
  guess <- tolower(paste0(substr(sname[1], 1, 1), sname[2]))
  if(paste0(guess, '_gene_ensembl') %in% dataset$dataset){
    return(switch(output,
                  'abbr'=guess,
                  'mart'=useEnsembl(biomart = 'ensembl',
                                    dataset=paste0(guess, '_gene_ensembl')),
                  'scientific name'=sname0,
                  'taxid'=taxonomy[id[1], 'taxid'],
                  'common name'=taxonomy[id[1], 'common_name']))
  }else{
    warning('The dataset is not available for ', sname)
  }
}
