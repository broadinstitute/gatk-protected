source("http://bioconductor.org/biocLite.R")
biocLite("DNAcopy")
argparserurl="http://cran.r-project.org/src/contrib/Archive/argparser/argparser_0.1.tar.gz"
if (!("argparser" %in% rownames(installed.packages))) {
  install.packages(argparserurl, repos=NULL, type="source")
}
dependencies = c("naturalsort")
if (!all(dependencies %in% rownames(installed.packages))) {
  install.packages(setdiff(dependencies, rownames(installed.packages())), repos="http://cran.cnr.Berkeley.edu")
}
q(save="no")
