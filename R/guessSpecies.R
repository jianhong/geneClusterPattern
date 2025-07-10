#' Guess the species
#' @description
#' Guess the species for a given string.
#' @param species A character(1L) to guess.
#' @param output Output type.
#' @param ... Parameter will be passed to \link[biomaRt]{useEnsembl}.
#' @return A character or Mart object.
#' @importFrom biomaRt useEnsembl listDatasets
#' @importFrom utils adist
#' @export
#' @examples
#' if(interactive()){
#'   guessSpecies('zebrafish')
#'   guessSpecies('Japanese medaka')
#'   guessSpecies('stickleback')
#'   guessSpecies('hsapiens', output='scientific name')
#'   guessSpecies('human', output='taxid')
#'   guessSpecies('hg38', output='scientific name')
#'   guessSpecies('hg38', output='taxid')
#'   guessSpecies('GRCh38', output='scientific name')
#'   guessSpecies('GRCh38', output='taxid')
#' }
#'
guessSpecies <- function(species,
                         output=c('abbr', 'scientific name',
                                  'mart',
                                  'taxid', 'common name'),
                         ...){
  output <- match.arg(output)
  if(length(species)>1){
    return(unlist(lapply(species, guessSpecies, output=output, ...)))
  }
  stopifnot(is.character(species))
  # load taxonomy from sysdata.rda
  if(paste0(species, '_gene_ensembl') %in% ensembl_dataset$dataset){
    id <- which(taxonomy$abbr==species & !is.na(taxonomy$abbr))
    return(switch(output,
                  'abbr'=species,
                  'mart'=useEnsembl(biomart = 'ensembl',
                                    dataset=paste0(species, '_gene_ensembl'),
                                    ...),
                  'scientific name'=taxonomy[id[1], 'scientific_name'],
                  'taxid'=taxonomy[id[1], 'tax_id'],
                  'common name'=taxonomy[id[1], 'common_name']))
  }
  if(species %in% ucsc_release$UCSC.VERSION){
    id <- which(ucsc_release$UCSC.VERSION==species & 
                  !is.na(ucsc_release$UCSC.VERSION))
    species <- ucsc_release$abbr[id]
    return(switch(output,
                  'abbr'=species[1],
                  'mart'=useEnsembl(biomart = 'ensembl',
                                    dataset=paste0(species, '_gene_ensembl'),
                                    ...),
                  'scientific name'=ucsc_release[id[1], 'scientific_name'],
                  'taxid'=ucsc_release[id[1], 'tax_id'],
                  'common name'=ucsc_release[id[1], 'common_name']))
  }
  if(species %in% ensembl_release$assembly){
    id <- which(ensembl_release$assembly==species & 
                  !is.na(ensembl_release$assembly))
    species <- ensembl_release$abbr[id]
    return(switch(output,
                  'abbr'=species[1],
                  'mart'=useEnsembl(biomart = 'ensembl',
                                    dataset=paste0(species, '_gene_ensembl'),
                                    ...),
                  'scientific name'=ensembl_release[id[1], 'scientific_name'],
                  'taxid'=ensembl_release[id[1], 'tax_id'],
                  'common name'=ensembl_release[id[1], 'common_name']))
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
  }else if(abs(min_scientific_name - min_common_name)<2){
    id <- which(dist_scientific_name==min_scientific_name |
                  grepl(species, taxonomy$scientific_name))
    id1 <- which(dist_common_name==min_common_name  |
                  grepl(species, taxonomy$common_name))
    out <- pasteReplaceLast(c(taxonomy[id, 'scientific_name'],
                              taxonomy[id1, 'common_name']))
    stop('Can not find exactly match for ', species,
         '. Do you mean ', out, '?\n')
  }else if(min_scientific_name < min_common_name){
    id <- which(dist_scientific_name==min_scientific_name |
                  grepl(species, taxonomy$scientific_name))
    out <- pasteReplaceLast(taxonomy[id, 'scientific_name'])
    stop('Can not find exactly match for ', species,
         '. Do you mean ', out, '?\n')
  }else{
    id <- which(dist_common_name==min_common_name  |
                   grepl(species, taxonomy$common_name))
    out <- pasteReplaceLast(taxonomy[id, 'common_name'])
    stop('Can not find exactly match for ', species,
         '. Do you mean ', out, '?\n')
  }
  sname0 <- taxonomy[id[1], 'scientific_name']
  sname <- strsplit(sname0, ' ')[[1]]
  guess <- tolower(paste0(substr(sname[1], 1, 1), sname[2]))
  if(paste0(guess, '_gene_ensembl') %in% ensembl_dataset$dataset){
    return(switch(output,
                  'abbr'=guess,
                  'mart'=useEnsembl(biomart = 'ensembl',
                                    dataset=paste0(guess, '_gene_ensembl'),
                                    ...),
                  'scientific name'=sname0,
                  'taxid'=taxonomy[id[1], 'tax_id'],
                  'common name'=taxonomy[id[1], 'common_name']))
  }else{
    stop('The dataset is not available for ', sname0)
  }
}
