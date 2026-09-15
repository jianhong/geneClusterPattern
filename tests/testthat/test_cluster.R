# prepare the test data
fish <- readRDS(system.file('extdata', 'fish.rds',
                            package = 'geneClusterPattern'))
homologs <- readRDS(system.file('extdata', 'homologs.rds',
                                package = 'geneClusterPattern'))
queryGene <- 'inhbaa'
queryGeneEnID <- 'ENSDARG00000012671'
targetGeneEnID <- c('ENSDARG00000102341')
queryGene2 <- 'pcolce2b'
queryGeneEnID2 <- 'ENSDARG00000055575'
targetGeneEnID2 <- c("ENSDARG00000055575", "ENSDARG00000061203",
                     "ENSDARG00000011821", "ENSDARG00000077468",
                     "ENSDARG00000031307", "ENSDARG00000015567",
                     "ENSDARG00000075369", "ENSDARG00000045447")
test_that("getGeneCluster works not correct", {
  nearestNeighbors <- getGeneCluster(fish, queryGene, homologs, k=3)
  expect_true(queryGeneEnID %in% nearestNeighbors)
  expect_true(targetGeneEnID %in% nearestNeighbors)
  expect_true(length(nearestNeighbors)==3)
  
  nearestNeighbors <- getGeneCluster(fish, queryGene2, homologs, k=10)
  expect_true(queryGeneEnID2 %in% nearestNeighbors)
  expect_all_true(targetGeneEnID2 %in% nearestNeighbors)
  expect_true(length(nearestNeighbors)==10)
})

test_that("getGeneRank works not correct", {
  ggr <- getGeneRank(fish)
  expect_all_true(names(ggr) %in% seqlevels(fish))
  expect_true(sum(lengths(ggr))==length(fish))
  chr1 <- fish[seqnames(fish)=="1"]
  strand(chr1) <- '*'
  chr1 <- chr1[order(start(chr1))]
  expect_identical(names(chr1), names(ggr[['1']]))
})

test_that("guessSpecies not work correct", {
  output <- c('abbr', 'scientific name',
              'taxid', 'common name')
  x <- unlist(lapply(output, guessSpecies, species='zebrafish'))
  expect_equal(c('drerio', 'Danio rerio', '7955', 'zebrafish'),
                   x)
  expect_equal(tolower(guessSpecies('hg38', output = 'common name')),
               'human')
  expect_equal(tolower(guessSpecies('GRCh38', output = 'common name')),
               'human')
  expect_equal(tolower(guessSpecies('mm10', output = 'common name')),
               'mouse')
  expect_equal(tolower(guessSpecies('GRCm39', output = 'common name')),
               'mouse')
})

test_that("homologsFromEnsemblIDs input check", {
  expect_error(homologsFromEnsemblIDs(mart='abc', ensembl_gene_ids='abc'))
  expect_error(homologsFromEnsemblIDs(ensembl_gene_ids=1))
})

test_that("grangesFromEnsemblIDs input check", {
  expect_error(grangesFromEnsemblIDs(mart='abc', ensembl_gene_ids='abc'))
  expect_error(grangesFromEnsemblIDs(ensembl_gene_ids=1))
})

test_that('addGeneInfo',{
  gr <- fish[targetGeneEnID2]
  id <- seq.int(7)
  colors <- setNames(id, targetGeneEnID2[id])
  x <- addGeneInfo(gr, colors = colors)
  expect_equal(x$color[id], palette()[id])
  expect_all_true(x$must_have_label[id])
  expect_all_equal(x$color[-id], '#33333311')
  expect_all_false(x$hide_label[id])
  expect_true(length(x$feature)==length(x))
})

test_that('getGeneClusterPattern', {
  gr <- fish[targetGeneEnID2]
  id <- seq.int(7)
  colors <- setNames(id, targetGeneEnID2[id])
  gr <- addGeneInfo(gr, colors = colors)
  reg <- GRanges('24:3347510-5910811')
  gr1 <- getGeneClusterPattern(gr, reg=reg, k=3)
  gr1 <- gr1[names(gr)]
  for(i in names(mcols(gr))){
    if(i!='feature'){
      expect_equal(mcols(gr)[[i]], mcols(gr1)[[i]])
    }
  }
  expect_all_false(start(gr1)==start(gr))
})

test_that('rescaleRegion', {
  gr <- fish[targetGeneEnID2]
  x <- rescaleRegion(gr)
  expect_all_equal(round(width(x)/max(width(gr)), digits = 3), 1)
})

test_that('getGeneModel', {
  gr <- fish[targetGeneEnID2]
  gr <- subsetByOverlaps(fish, range(gr))
  ggm <- getGeneModel(gr, ids=targetGeneEnID2, max_gap = 1e7)
  expect_all_true(c('geneModel', 'region', 'id0_width_pct', 'id_width_pct') %in%
                    names(ggm))
  expect_true(ggm$id0_width_pct < ggm$id_width_pct)
  id0_width_pct <- width(gr[targetGeneEnID2[1]])/width(ggm$region)
  expect_equal(id0_width_pct, ggm$id0_width_pct)
})

test_that('checkAndGetGeneModels', {
  gr <- fish[targetGeneEnID2]
  genesList <- GRangesList(gr)
  expect_error(checkAndGetGeneModels(genesList, 
                                     ids = targetGeneEnID2,
                                     max_gap = 1e7),
               "genesList must have names")
  genesList <- GRangesList(A=gr)
  genesList <- lapply(genesList, function(.ele){
    .ele$gene_name <- NULL
    .ele
  })
  expect_error(checkAndGetGeneModels(genesList, 
                                     ids = targetGeneEnID2,
                                     max_gap = 1e7),
               'Elements in genesList must have metadata "gene_name"')
  genesList <- GRangesList(A=gr)
  genesList <- lapply(genesList, unname)
  expect_error(checkAndGetGeneModels(genesList, 
                                     ids = targetGeneEnID2,
                                     max_gap = 1e7),
               'Elements in genesList must have names')
})

test_that('pruningSequences', {
  # remove chrom from 10 to 25
  x <- pruningSequences(list(fish), negativePattern = '\\d{2}')
  expect_all_true(seqlevels(x[[1]]) %in% as.character(seq.int(9)))
})

test_that('orthologPairsFromOrthoFinder', {
  path <- system.file('extdata/orthofinder/Danio_rerio.GRCz11.pep.all.tsv.gz',
                    package='geneClusterPattern')
  orthologs <- orthologPairsFromOrthoFinder(path)
  expect_true(length(orthologs)==5)
  n <- vapply(orthologs, ncol, numeric(1L))
  expect_all_equal(n, 2)
})

test_that('getHomologListForOrthoFinder', {
  path <- system.file('extdata/orthofinder/Danio_rerio.GRCz11.pep.all.tsv.gz',
                      package='geneClusterPattern')
  orthologs <- orthologPairsFromOrthoFinder(path)
  ## annotations
  annoGR_list <- readRDS(
    system.file('extdata/orthofinder/grange.obj.rds',
                package='geneClusterPattern'))
  ## reformat to meet the requirement of the package
  homologs <- getHomologListForOrthoFinder(orthologs, annoGR_list[-1])
  n <- lapply(homologs, function(.ele){
    expect_true(length(names(.ele))==length(.ele))
    expect_all_true(c('gene_name', 'homolog_ensembl_gene_ids') %in%
                      colnames(mcols(.ele)))
    expect_all_true(grepl('^ENSDARG', names(.ele)))
    expect_all_true(grepl('^ENS', .ele$homolog_ensembl_gene_ids))
  })
})


test_that('alignCenterId', {
  geneClusterPatterns <- list(
    species1 = GRanges('chr1', IRanges(start = c(10, 30, 50),
                                       end = c(20, 40, 60),
                                       names = c("geneA", "geneB", "geneC"))),
    species2 = GRanges('XII', IRanges(start = c(5, 25, 70),
                                      end = c(20, 35, 80),
                                      names = c("geneA", "geneC", "geneB")))
  )
  patternRegions <- lapply(geneClusterPatterns, function(gr){
    strand(gr) <- '*'
    reg <- range(gr)
    ranges(reg) <- IRanges(start = 1, end = end(reg)+3)
    return(reg)
  })
  ## make the plot region size comparable
  patternRegions <- rescaleRegion(patternRegions)
  ## by local
  res <- alignCenterId(geneClusterPatterns, patternRegions, id0 = "geneA")
  ## test the output 
  expect_type(res, "list")
  expect_named(res, c("patterns", "regions"))
  # geneA center in species1
  ctr <- function(x) start(x) + width(x)/2
  # geneA should now sit exactly at the center of its (widened) region
  for (i in seq_along(res$patterns)) {
      pattern <- res$patterns[[i]]
      center_gene   <- ctr(pattern['geneA'])
      center_region <- width(res$regions)[i] / 2
      expect_equal(center_gene, center_region, tolerance = 0.5)
      # relative position of geneB should not change
      expect_equal(start(pattern['geneB']),
                   end(pattern['geneA'])+
                     distance(geneClusterPatterns[[i]]['geneB'],
                              geneClusterPatterns[[i]]['geneA'])+1)
  }
  ## global
  res <- alignCenterId(geneClusterPatterns, patternRegions)
  expect_type(res, "list")
  expect_named(res, c("patterns", "regions"))
  expect_equal(ctr(unlist(range(GRangesList(res$patterns)))),
               ctr(res$regions),
               tolerance = 0.5)
  for (i in seq_along(res$patterns)) {
    pattern <- res$patterns[[i]]
    # relative position of geneB should not change
    expect_equal(start(pattern['geneB']),
                 end(pattern['geneA'])+
                   distance(geneClusterPatterns[[i]]['geneB'],
                            geneClusterPatterns[[i]]['geneA'])+1)
  }
})

test_that('geneOrderScore', {
  genesList <- list(
    species1 = GRanges('chr1', IRanges(start = c(10, 30, 50),
                                       end = c(20, 40, 60),
                                       names = c("geneA", "geneB", "geneC")),
                       strand = c('+', '+', '-'),
                       gene_name=c("geneA", "geneB", "geneC")),
    species2 = GRanges('XII', IRanges(start = c(5, 25, 70),
                                      end = c(20, 35, 80),
                                      names = c("geneA", "geneC", "geneB")),
                       strand = c('+', '-', '+'),
                       gene_name=c("geneA", "geneC", "geneB"))
  )
  ids <- c('geneA', 'geneB', 'geneC')
  gos <- geneOrderScore(genesList, ids, method='edit distance')
  expect_equal(gos, 1/(adist('ABC', 'ACB')[1]+1))
  gos <- geneOrderScore(genesList, ids, method='edit distance', 
                        output = "value matrix")
  expect_equal(gos, matrix(c(1, 1/3, 1/3, 1), nrow = 2,
                           dimnames = list(names(genesList),
                                           names(genesList))))

  gos <- geneOrderScore(genesList, ids, method='global alignment score')
  expect_equal(gos, pwalign::score(pairwiseAlignment('ABC', 'ACB')))
  
  gos <- geneOrderScore(genesList, ids, method='spearman correlation')
  expect_equal(gos, cor(c(1, 2, 3), c(1, 3, 2), method = 'spearman'))
  
  gos <- geneOrderScore(genesList, ids, method='kendall correlation')
  expect_equal(gos, cor(c(1, 2, 3), c(1, 3, 2), method = 'kendall'))
  
  gos <- geneOrderScore(genesList, ids, method='non-random score')
  expect_equal(gos, jaccard(
    c('A', 'B', 'C', 'AB', 'BC', 'ABC'),
    c('A', 'C', 'B', 'AC', 'BC', 'ABC')
  ))
  
  gos <- geneOrderScore(genesList, ids, method='pairs direction')
  expect_equal(gos, 1) ## all strand pairs are identical
  
  strand(genesList$species1)[3] <- '+'
  gos <- geneOrderScore(genesList, ids, method='pairs direction')
  expect_equal(gos, 1/3) ## only AB is same stranded, AC and BC not same
})

test_that('pasteReplaceLast', {
  expect_equal(pasteReplaceLast(character(0)), "")
  expect_equal(pasteReplaceLast("apple"), "apple")
  expect_equal(pasteReplaceLast(c("apple", "banana")), "apple, or banana")
  expect_equal(
    pasteReplaceLast(c("apple", "banana", "cherry")),
    "apple, banana, or cherry"
  )
  expect_equal(
    pasteReplaceLast(c("apple", "banana", "cherry"), collapse = "; "),
    "apple; banana; cherry"
  )
  expect_equal(
    pasteReplaceLast(c("apple", "banana", "cherry"), last = ", and"),
    "apple, banana, and cherry"
  )
})

test_that('convertNamesOfHomologIDs', {
  homologs <- list(
    species1 = GRanges('chr1', IRanges(start = c(10, 30, 50),
                                       end = c(20, 40, 60),
                                       names = c("geneA", "geneB", "geneC")),
                       strand = c('+', '+', '-'),
                       homolog_ensembl_gene_ids=c("A", "B", "C")),
    species2 = GRanges('XII', IRanges(start = c(5, 25, 70),
                                      end = c(20, 35, 80),
                                      names = c("geneA", "geneC", "geneB")),
                       strand = c('+', '-', '+'),
                       homolog_ensembl_gene_ids=c("A", "C", "B")),
    species3 = GRanges()
  )
  hl <- convertNamesOfHomologIDs(homologs)
  for(i in seq_along(homologs)){
    expect_equal(names(hl[[i]]), homologs[[i]]$homolog_ensembl_gene_ids)
    expect_equal(hl[[i]]$homolog_ensembl_gene_ids, names(homologs[[i]]))
  }
})
