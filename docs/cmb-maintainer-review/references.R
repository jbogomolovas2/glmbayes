# Reuse unchanged numerical integration code only; no sampler run here.
library(glmbayes)
args<-commandArgs(TRUE);name<-args[1];out<-'references';dir.create(out,showWarnings=FALSE)
lines<-readLines('source/inst/validation/cmb-matrix.R')
start<-grep('^set.seed',lines)[1];end<-grep('^rows<-list',lines)[1]-1
print(sessionInfo());cat('mvtnorm',as.character(packageVersion('mvtnorm')),'\n');eval(parse(text=lines[start:end]))
