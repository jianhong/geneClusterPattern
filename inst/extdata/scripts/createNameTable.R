download.file('https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz', 'taxdump.tar.gz')
untar('taxdump.tar.gz')
names <- read.delim('names.dmp', header=FALSE, sep = '\t')
names <- names[, c(1, 3, 5, 7)]
colnames(names) <- c('tax_id', 'name_txt', 'unique_name', 'name_class')
names[grepl('human', names[, 2]), ]

division <- read.delim('division.dmp', header=FALSE, sep='\t')
divisions <- division[!division$V5 %in% c('Bacteria', 'Phages', 'Synthetic and Chimeric', 'Unassigned', 'Viruses', 'Environmental samples'), 1]

nodes <- read.delim('nodes.dmp', header=FALSE, sep = '\t')
nodes <- nodes[, c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25)]
colnames(nodes) <- c('tax_id', 'parent tax_id', 'rank', 'embl code',
                     'division id', 'inherited div flag',
                     'genetic code id', 'inherited GC  flag',
                     'mitochondrial genetic code id',
                     'inherited MGC flag',
                     'GenBank hidden flag',
                     'hidden subtree root flag',
                     'comments')
nodes <- nodes[, c('tax_id', 'parent tax_id', 'division id')]
nodes <- unique(nodes)

tax_id <- nodes[nodes[, 'division id'] %in% divisions, 'tax_id']
tax_id <- unique(tax_id)

keep <- names[names$tax_id %in% tax_id, ]
keep <- keep[keep$name_class %in% c('scientific name', 'genbank common name'), ]
keep <- keep[!grepl("[^a-zA-Z ]", keep[, 2]), ]
dim(keep)
keep[grep('zebrafish', keep[, 2]), ]
tax_id_len <- table(keep$tax_id)
keep <- keep[keep$tax_id %in% names(tax_id_len[tax_id_len==2]), ]
keep_s <- keep[keep$name_class=='scientific name', c(1, 2)]
keep_c <- keep[keep$name_class=='genbank common name', c(1, 2)]
colnames(keep_s)[2] <- 'scientific_name'
colnames(keep_c)[2] <- 'common_name'
name <- merge(keep_s, keep_c)
name <- name[grepl(' ', name$scientific_name), ]
saveRDS(name, 'taxonomy_names.rds')
abbr <- name$scientific_name
abbr <- strsplit(abbr, '\\s+')
keep <- lengths(abbr)==2
name <- name[keep, ]
abbr <- name$scientific_name
abbr <- strsplit(abbr, '\\s+')
abbr <- do.call(rbind, abbr)
name$abbr <- tolower(paste0(substr(abbr[, 1], 1, 1), abbr[, 2]))
taxonomy <- name
mart <- useEnsembl(biomart = 'ensembl')
ensembl_dataset <- listDatasets(mart = mart)
use_data(taxonomy, ensembl_dataset, internal = TRUE, overwrite = TRUE, compress = 'xz')

unlink(list.files(pattern='dmp'))
unlink('gc.prt')
unlink('readme.txt')
unlink('taxdump.tar.gz')

library(geneClusterPattern)
library(org.Dr.eg.db)
## prepare all the ensembl ids
ids <- as.list(org.Dr.egENSEMBL)
ensembl_gene_ids <- sort(unique(unlist(ids)))
fish_mart <- guessSpecies('zebrafish', output='mart')
fish <- grangesFromEnsemblIDs(mart = fish_mart,
                              ensembl_gene_ids = ensembl_gene_ids)
ensembl_gene_ids <- names(fish[seqnames(fish)=='24'])
## keep the standard sequence only
fish <- pruningSequences(fish)
saveRDS(fish, 'inst/extdata/fish.rds')

species <- guessSpecies(c('human', 'house mouse', 'Japanese medaka', 'turquoise killifish', 'gaculeatus')) # three-spined stickleback
homologs <- getHomologGeneList(species, fish_mart, ensembl_gene_ids)
homologs <- pruningSequences(homologs)
saveRDS(homologs, 'inst/extdata/homologs.rds')

# ucsc_release <- read.csv('ignore/ucsc.release.versions.2025.07.csv')
library(httr)
library(jsonlite)
library(xml2)
for(server in c('may2025', 'oct2024', 'may2024', 'jan2024', 'jul2023',
                'feb2023', 'oct2022', 'jul2022', 'apr2022', 'dec2021',
                'may2021', 'feb2021', 'nov2020', 'aug2020', 'apr2020',
                'jan2020', 'may2015', 'oct2014', 'feb2014', 'may2009')){
  r <- GET(paste0("https://", server, ".rest.ensembl.org/info/species?"),
           content_type("application/json"))
  if(r$status_code==200){
    r <- fromJSON(toJSON(content(r)))[[1]]
    r1 <- r[c("assembly", "common_name", "release", "display_name", "name", "taxon_id")]
    r1 <- lapply(r1, function(.ele){
      sapply(.ele, `[`, i=1)
    })
    r1 <- do.call(cbind, r1)
    r1 <- as.data.frame(r1)
    write.csv(r1, paste0('ignore/ensembl.',server,'.csv'), row.names=FALSE)
  }
}
for(server in c('may2025', 'oct2024', 'may2024', 'jan2024', 'jul2023',
                'feb2023', 'oct2022', 'jul2022', 'apr2022', 'dec2021',
                'may2021', 'feb2021', 'nov2020', 'aug2020', 'apr2020',
                'jan2020', 'may2015', 'oct2014', 'feb2014', 'may2009')){
  r <- GET(paste0("https://", server, ".rest.ensembl.org/info/species?"),
           content_type("application/json"))
  if(r$status_code==200){
    r <- fromJSON(toJSON(content(r)))[[1]]
    r1 <- r[c("assembly", "common_name", "release", "display_name", "name", "taxon_id")]
    r1 <- lapply(r1, function(.ele){
      sapply(.ele, `[`, i=1)
    })
    r1 <- do.call(cbind, r1)
    r1 <- as.data.frame(r1)
    write.csv(r1, paste0('ignore/ensembl.',server,'.csv'), row.names=FALSE)
  }
}

ensembl_release <- lapply(dir('ignore', 'ensembl', full.names = TRUE), read.csv)
ensembl_release <- do.call(rbind, ensembl_release)
write.csv(ensembl_release, 'ignore/ensembl_release.version.202507.csv', row.names = FALSE)

r <-  GET('https://api.genome.ucsc.edu/list/ucscGenomes')
if(r$status_code==200){
  r <- fromJSON(toJSON(content(r)))[['ucscGenomes']]
  n <- unique(unlist(lapply(r, names)))
  r <- do.call(rbind, lapply(r, function(.ele) .ele[n]))
  write.table(r, 'ignore/ucsc_release.versions.2025.07.tab', row.names = TRUE, quote = TRUE, sep = '\t', col.names = TRUE)
}

rm(list=ls())
load('R/sysdata.rda')
ensembl_release <- read.csv('ignore/ensembl_release.version.202507.csv')
abbr <- do.call(rbind, strsplit(tolower(ensembl_release$name), '_'))
abbr <- paste0(substr(abbr[, 1], 1, 1), abbr[, 2])
ensembl_release$abbr <- abbr
colnames(ensembl_release)[5:6] <- c('scientific_name', 'tax_id')

ucsc_release <- read.delim('ignore/ucsc_release.versions.2025.07.tab')
ucsc_release$UCSC.VERSION <- rownames(ucsc_release)
abbr <- do.call(rbind, strsplit(tolower(ucsc_release$scientificName), ' '))
abbr <- paste0(substr(abbr[, 1], 1, 1), abbr[, 2])
ucsc_release$abbr <- abbr
ucsc_release <- ucsc_release[, c('taxId', 'scientificName', 'organism', 'abbr', 'UCSC.VERSION', 'description', 'genome')]
colnames(ucsc_release)[1:4] <- c('tax_id', 'scientific_name','common_name','abbr')
rm(abbr)

save.image('R/sysdata.rda', compress='xz')
