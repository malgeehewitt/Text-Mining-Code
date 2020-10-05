getVariableWords<-function(dtm, threshold=1){
  std.devs<-NULL
  words<-colnames(dtm)
  for (i in 1:ncol(dtm)){
    curr.stdev<-sd(dtm[,i])
    std.devs<-c(std.devs, curr.stdev)
  }
  mean.sd<-mean(std.devs)
  sd.sd<-sd(std.devs)
  
}

window<-function(corpus.text, window.size, window.offset){
  window.vector<-NULL
  text.words<-unlist(strsplit(corpus.text, " "))
  max.window.start<-length(text.words)-window.size
  num.windows<-max.window.start %/% window.offset
  window.start<-1
  window.end<-window.start+window.size
  for (i in 1:num.windows){
    if (i == num.windows){
      curr.window<-text.words[window.start:length(text.words)]
    } else {
      curr.window<-text.words[window.start:window.end]
      window.start<-window.start+window.offset
      window.end<-window.start+window.size
    }
    curr.window<-paste(curr.window, collapse=" ")
    window.vector<-c(window.vector, curr.window)
  }
  window.vector<-VectorSource(window.vector)
  windowed.corpus<-Corpus(window.vector, readerControl=list(language="English"))
}

findIncDecWords<-function(windowed.dtm, threshold=2){
  all.slopes<-NULL
  inc.slopes<-NULL
  dec.slopes<-NULL
  xvalues<-seq(1:nrow(windowed.dtm))
  words<-colnames(windowed.dtm)
  for (i in 1:ncol(windowed.dtm)){
    curr.word<-windowed.dtm[,i]
    mydata<-cbind(xvalues, curr.word)
    mydata<-as.data.frame(mydata)
    colnames(mydata)<-c("x","y")
    curr.word.lm<-lm(y~x, data=mydata)
    curr.slope<-coef(curr.word.lm)[2]
    all.slopes<-c(all.slopes, curr.slope)
  }
  names(all.slopes)<-words
  slope.sd<-sd(all.slopes)
  slope.mean<-mean(all.slopes)
  inc.slope.index<-which(all.slopes>(slope.mean+(threshold*slope.sd)))
  dec.slope.index<-which(all.slopes<(slope.mean-(threshold*slope.sd)))
  #print (inc.slope.index)
  #print (dec.slope.index)
  inc.slopes<-all.slopes[inc.slope.index]
  dec.slopes<-all.slopes[dec.slope.index]
  sig.slopes<-list(inc.slopes, dec.slopes)
  names(sig.slopes)<-c("increasing", "decreasing")
  return(sig.slopes)
}

findVWords<-function(windowed.dtm, threshold=2){
  all.sds<-NULL
  v.words<-NULL
  xvalues<-seq(1:nrow(windowed.dtm))
  words<-colnames(windowed.dtm)
  for (i in 1:ncol(windowed.dtm)){
    curr.word<-windowed.dtm[,i]
    curr.word.sd<-sd(curr.word)
    all.sds<-c(all.sds, curr.word.sd)
  }
  names(all.sds)<-words
  sd.sds<-sd(all.sds)
  sd.mean<-mean(all.sds)
  v.index<-which(all.sds>(sd.mean+(threshold*sd.sds)))
  #print (inc.slope.index)
  #print (dec.slope.index)
  v.words<-all.sds[v.index]
  v.words<-sort(v.words, decreasing=T)
  return(v.words)
}
