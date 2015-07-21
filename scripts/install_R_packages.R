source("http://bioconductor.org/biocLite.R")
biocLite("DNAcopy")
optparseurl="http://cran.r-project.org/src/contrib/optparse_1.3.0.tar.gz"
if (!("optparse" %in% rownames(installed.packages))) {
  install.packages(optparseurl, repos=NULL, type="source")
}

dependencies = c("naturalsort")
if (!all(dependencies %in% rownames(installed.packages))) {
  install.packages(setdiff(dependencies, rownames(installed.packages())), repos="http://cran.cnr.Berkeley.edu")
}
q(save="no")
