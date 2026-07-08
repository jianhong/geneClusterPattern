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
#            Danio rerio Gasterosteus aculeatus   Lepisosteus oculatus 
#                  37241                  30416                  23315 
# Nothobranchius furzeri        Oryzias latipes      Takifugu rubripes 
#                  25475                  24365                  24406 
bm <- GRangesList(bm)

## only keep required chromosomes 24 for fish
bm_zfish <- bm$`Danio rerio`[seqnames(bm$`Danio rerio`) %in% '24']
fish2all <- read.delim('ignore/orthofinder/Danio_rerio.GRCz11.pep.all.tsv.gz')
keep <- lapply(names(bm_zfish), function(id) grepl(id, fish2all$Danio_rerio.GRCz11.pep.all))
keep <- do.call(cbind, keep)
keep <- rowSums(keep)>0
table(keep)
fish2all.sub <- fish2all[keep, ]
write.table(fish2all.sub, gzfile('inst/extdata/orthofinder/Danio_rerio.GRCz11.pep.all.tsv.gz'), quote = FALSE, sep = '\t', row.names = FALSE)

bm_to_keep <- with(fish2all.sub, split(Orthologs, Species))
bm_to_keep <- lapply(bm_to_keep, function(.ele) 
  trimENSname(unlist(strsplit(.ele, ', '))))
names(bm_to_keep)
names(bm)
bm.sub <- bm
bm.sub[[1]] <- bm_zfish
for(i in seq_along(bm)[-1]){
  bm.sub[[i]] <- bm[[i]][names(bm[[i]]) %in% bm_to_keep[[i-1]]]
}
lengths(bm.sub)
#            Danio rerio Gasterosteus aculeatus   Lepisosteus oculatus 
#                    849                    590                    571 
# Nothobranchius furzeri        Oryzias latipes      Takifugu rubripes 
#                    525                    603                    568 
saveRDS(bm.sub, 'inst/extdata/orthofinder/grange.obj.rds')

orthogroups <- read.delim('ignore/orthofinder/Orthogroups.tsv.gz')
keep <- lapply(names(bm_zfish), function(id) grepl(id, orthogroups$`Danio rerio`))
keep <- do.call(cbind, keep)
keep <- rowSums(keep)>0
table(keep)
orthogroups.sub <- orthogroups[keep, ]
write.table(orthogroups.sub, gzfile('inst/extdata/orthofinder/Orthogroups.tsv.gz'),
            quote = FALSE, sep = '\t', row.names = TRUE)

