# Don't activate renv in an OpenScPCA docker image
if (Sys.getenv('OPENSCPCA_DOCKER') != 'TRUE') {
  source('renv/activate.R')
}

# work around renv 1.0.11 bug
if (exists("object")){
  rm(object)
}
