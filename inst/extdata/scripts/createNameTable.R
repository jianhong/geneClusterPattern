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

unlink(list.files(pattern='dmp'))
unlink('gc.prt')
unlink('readme.txt')
unlink('taxdump.tar.gz')
