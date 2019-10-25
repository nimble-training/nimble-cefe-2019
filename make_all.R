setwd('~/github/nimble/nimble-cefe-2019')
library(methods)
library(rmarkdown)
##
(allModules <- gsub('\\.Rmd$', '', list.files('docs', pattern = '*\\.Rmd')))
dirname <- 'docs'
##
for(module in allModules) {
    message('processing: ', module)
    rmdFile <- paste0(dirname, '/', module, '.Rmd')
    htmlFile <- paste0(dirname, '/', module, '.html')
    system(paste0('rm -rf ', htmlFile), )
    generatedFile <- render(rmdFile, quiet = TRUE)
    if(!is.character(generatedFile)) stop('generated file name from render is not a character string, for module: ', module)
    if(!file.exists(generatedFile)) stop('output html file doesn\'t exist, for module: ', module)
}





