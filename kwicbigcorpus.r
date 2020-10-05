cleanText<-function(text.paste){
  text.split<-unlist(strsplit(text.paste, ''))
  text.split<-text.split[which(text.split %in% c(letters, LETTERS, ' '))]
  text.split<-tolower(text.split)
  text.paste<-paste(text.split, collapse='')
  text.paste<-unlist(strsplit(text.paste, ' '))
  return(text.paste)
}

getKwics<-function(text.split, target, term, horizon, context){
  kwic.hits<-which(text.split==target)
  if(length(kwic.hits)>0){
    start.pos<-kwic.hits-horizon
    end.pos<-kwic.hits+horizon
    bad.start<-which(start.pos<1)
    if(length(bad.start)>0){
      start.pos[bad.start]<-1
    }
    bad.end<-which(end.pos>length(text.split))
    if(length(bad.end)>1){
      end.pos[bad.end]<-length(text.split)
    }
    all.kwics<-mapply(function(x,y) text.split[x:y], start.pos, end.pos, SIMPLIFY = F)
    term.hits<-unlist(lapply(all.kwics, function(x) length(which(x==term))))
    term.hits<-which(term.hits>0)
    if(length(term.hits)>0){
      kwic.hits<-kwic.hits[term.hits]
      context.start<-kwic.hits-(horizon+context)
      context.end<-kwic.hits+(horizon+context)
      bad.con.start<-which(context.start<1)
      if(length(bad.con.start)>0){
        context.start[bad.con.start]<-1
      }
      bad.con.end<-which(context.end>length(text.split))
      if(length(bad.con.end)>0){
        context.end[bad.con.end]<-length(text.split)
      }
      context.kwics<-unlist(mapply(function(x,y) paste(text.split[x:y], collapse=" "), context.start, context.end, SIMPLIFY=T))
      target.col<-rep(target, length(context.kwics))
      term.col<-rep(term, length(context.kwics))
      kwic.df<-data.frame(target.col, term.col, context.kwics)
      colnames(kwic.df)<-c("Target", "Term", "Kwics")
      return(kwic.df)
    }
  }
}

getTextKwics<-function(textname, term.list, target.list, horizon, context){
  print(textname)
  text.import<-scan(textname, what='character', sep="\n", quiet=T)
  text.import<-paste(text.import, collapse=" ")
  text.clean<-cleanText(text.import)
  kwic.tables<-mapply(function(x,y) getKwics(text.clean, x, y, horizon, context), term.list, target.list, SIMPLIFY = F)
  kwic.tables<-do.call('rbind', kwic.tables)
  return(kwic.tables)
}

kwicListPairs<-function(folder, target.list, term.list, horizon, context){
  if(length(target.list) != length(term.list)){
    print("Length of target and term vectors do not match")
    return(NULL)
  }
  all.files<-list.files(folder, pattern='.txt', full.names=T)
  all.kwic.tables<-lapply(all.files, function(x) getTextKwics(x, target.list, term.list, horizon, context))
  all.kwic.tables<-do.call('rbind', all.kwic.tables)
  #return(all.kwic.tables)
  all.kwic.tables<-all.kwic.tables[order(all.kwic.tables$Term),]
  all.kwic.tables<-all.kwic.tables[order(all.kwic.tables$Target),]
  write.csv(all.kwic.tables, file="KwicResults.csv", row.names=F)
}