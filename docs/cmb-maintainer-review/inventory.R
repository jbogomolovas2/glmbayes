# Parse function definitions without evaluating package code.
script_arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
review_dir<-dirname(normalizePath(sub('^--file=','',script_arg[[1]])))
base<-'7959a42d32f822a4b65ddfaf94927190ad5db536';head<-'b971596511c90f33de290497b2c7e8feb1a4e580'
files<-system2('git',c('diff','--name-only',base,head,'--','R'),stdout=TRUE)
collect<-function(expr,prefix='') {
 out<-list()
 walk<-function(e,path='') {
  if(!is.call(e))return()
  if(as.character(e[[1]])[1] %in% c('<-','=') && length(e)==3 && is.symbol(e[[2]]) && is.call(e[[3]]) && identical(e[[3]][[1]],as.name('function'))) {
   nm<-paste0(path,as.character(e[[2]]));out[[nm]]<<-e[[3]];walk(e[[3]][[3]],paste0(nm,'::'));return()
  }
  if(length(e)>1)for(i in seq.int(2,length(e))) {if(identical(e[[i]],quote(expr=)))next;if(!is.symbol(e[[i]]))walk(e[[i]],path)}
 }
 for(e in expr)walk(e,prefix)
 out
}
rows<-list()
for(f in files) {
 readrev<-function(rev) {s<-suppressWarnings(system2('git',c('show',paste0(rev,':',f)),stdout=TRUE,stderr=FALSE));if(!length(s))return(list());collect(parse(text=s))}
 old<-readrev(base);new<-readrev(head)
 for(nm in union(names(old),names(new))) {
  o<-old[[nm]];nn<-new[[nm]]
  if(identical(o,nn))next
  sig<-function(x)if(is.null(x))'absent' else paste(deparse(x[[2]],width.cutoff=500),collapse=' ')
  rows[[length(rows)+1]]<-data.frame(file=f,symbol=nm,old=sig(o),new=sig(nn),kind=if(is.null(o))'new' else if(identical(o[[2]],nn[[2]]))'body only' else 'signature')
 }
}
write.csv(do.call(rbind,rows),file.path(review_dir,'r-function-inventory.csv'),row.names=FALSE)
