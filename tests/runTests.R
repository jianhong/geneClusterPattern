require("geneClusterPattern") || stop("unable to load Package:geneClusterPattern")
require('S4Vectors') || stop('unable to load S4Vectors')
require("testthat") || stop("unable to load testthat")
test_check("geneClusterPattern")
