library(glmbayes)
# Reproduce the wrapper defect without changing the package.
err<-tryCatch(EnvelopeBuild(matrix(0,1,1),matrix(1,1,1),0,matrix(1,1,1),matrix(0,1,1),matrix(1,1,1),0,1,n_envopt=1),error=conditionMessage)
writeLines(as.character(err),'wrapper-defect.txt')
