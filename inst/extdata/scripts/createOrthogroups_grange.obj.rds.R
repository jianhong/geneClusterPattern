# check the mart
marts <- lapply(c("Danio rerio",
                  'Gasterosteus aculeatus',
                  'Lepisosteus oculatus',
                  'Nothobranchius furzeri',
                  'Oryzias latipes',
                  'Takifugu rubripes'),
                function(species) {
                    tryCatch({
                      guessSpecies(species, output='mart')
                    }, error=function(e){
                      message(e)
                      NULL
                    })
                  })
library(biomaRt)
mt <- listDatasets(mart=useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL"))
mt[grepl('stickleback', mt$description, ignore.case = TRUE), ]
GA_mt <- useMart('ensembl', 'gaculeatus_gene_ensembl')
names(marts) <- c("Danio rerio",
                  'Gasterosteus aculeatus',
                  'Lepisosteus oculatus',
                  'Nothobranchius furzeri',
                  'Oryzias latipes',
                  'Takifugu rubripes')
marts[[2]] <- GA_mt
# retrieve all gene positions
# bm <- lapply(marts, function(mart){ # biomart not stable
#   grangesFromEnsemblIDs(mart=mart)
# })
bm <- list()
for(i in seq_along(marts)){
  bm[[names(marts)[i]]] <- grangesFromEnsemblIDs(mart=marts[[i]])
}
lengths(bm)
     #       Danio rerio Gasterosteus aculeatus   Lepisosteus oculatus Nothobranchius furzeri        Oryzias latipes 
     #             37241                  30416                  23315                  25475                  24365 
     # Takifugu rubripes 
     #             24406 
bm <- GRangesList(bm)
saveRDS(bm, 'inst/extdata/orthofinder/grange.obj.rds')
