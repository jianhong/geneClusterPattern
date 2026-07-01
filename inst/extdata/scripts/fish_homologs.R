library(geneClusterPattern)
library(org.Dr.eg.db)
library(GenomeInfoDb)
ids <- as.list(org.Dr.egENSEMBL)
ensembl_gene_ids <- sort(unique(unlist(ids)))
## extract gene information via biomaRt
fish_mart <- guessSpecies('zebrafish', output='mart', version=112)
fish <- grangesFromEnsemblIDs(mart = fish_mart,
                              ensembl_gene_ids = ensembl_gene_ids)
saveRDS(fish, 'inst/extdata/fish.rds')

## subset the ensembl_gene_ids to save time
ensembl_gene_ids <- names(fish[seqnames(fish)=='24'])
## keep the standard sequence only
fish <- pruningSequences(fish)
## define the species scientific name
species <- guessSpecies(c('human', 'house mouse', 'Japanese medaka', 'turquoise killifish', 'gaculeatus'), version=112) # three-spined stickleback
species
homologs <- getHomologGeneList(species, fish_mart, ensembl_gene_ids)
homologs <- pruningSequences(homologs)
saveRDS(homologs, 'inst/extdata/homologs.rds')
