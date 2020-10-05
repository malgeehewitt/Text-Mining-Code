makeTopicModel<-function(text.corpus, ntopic){
  print(length(text.corpus))
  print("Making Topic Model")
  text.corpus<-tm_map(text.corpus, content_transformer(stemDocument))
  text.dtm<-DocumentTermMatrix(text.corpus, control=list(stopwords=T))
  #test.m<-as.matrix(text.dtm)
  #test.m<-rowSums(test.m)
  #print(test.m)
  text.topics<-LDA(text.dtm, ntopic, method="Gibbs")
  text.pp<-posterior(text.topics)
  text.pp<-text.pp$topics
  return(text.pp)
}

cleanText<-function(raw.text){
  clean.text<-unlist(strsplit(raw.text, ""))
  clean.text<-tolower(clean.text)
  remove.index<-which(clean.text %in% c(as.character(seq(0,9, by=1)), ",", ".", ";", ":", "'", '"', "?", "!", "&", "$"))
  clean.text<-clean.text[-remove.index]
  #keep.index<-which(clean.text %in% c(letters, " "))
  #clean.text<-clean.text[keep.index]
  clean.text<-paste(clean.text, collapse="")
  clean.text<-unlist(strsplit(clean.text, " "))
  return(clean.text)
}

divideTexts<-function(clean.text, div.length, text.name){
  print(text.name)
  text.starts<-seq(1, length(clean.text), by=div.length)
  if(length(text.starts)>10){
    text.ends<-text.starts[2:length(text.starts)]-1
    text.ends<-c(text.ends, length(clean.text))
    text.starts<-text.starts[1:(length(text.starts)-1)]
    text.ends<-text.ends[1:(length(text.ends)-1)]
    text.segs<-mapply(function(x,y) clean.text[x:y], text.starts, text.ends, SIMPLIFY=F)
    print(length(text.segs))
    text.segs<-lapply(text.segs, function(x) paste(x, collapse=" "))
    text.segs<-unlist(text.segs)
    return(text.segs)
  } else {
    return(NA)
  }
}

numTopic<-function(topic.row, cutoff=.5){
  topic.row<-sort(topic.row, decreasing=T)
  topic.sums<-unlist(lapply(seq(1,length(topic.row), by=1), function(x) sum(topic.row[1:x])))
  topic.overs<-which(topic.sums>cutoff)
  n.topic<-topic.overs[1]
  return(n.topic)
}

topicity<-function(clean.text, text.name){
  print(text.name)
  library(tm)
  library(topicmodels)
  text.corpus<-SimpleCorpus(VectorSource(clean.text))
  topic.pp<-makeTopicModel(text.corpus, 50)
  topicity<-NULL
  for(i in 1:nrow(topic.pp)){
    curr.topicity<-numTopic(topic.pp[i,])
    topicity<-c(topicity, curr.topicity)
  }
  seg.numbers<-as.character(seq(1,length(clean.text), by=1))
  seg.names<-rep(text.name, length(seg.numbers))
  topicity.table<-data.frame(seg.names, seg.numbers, topicity, stringsAsFactors=F)
  return(topicity.table)
}

sampleText<-function(text.segments, n.sample){
  if(length(text.segments)>n.sample){
    sampled.text<-sample(text.segments, n.sample)
    return(sampled.text)
  } else {
    return(NA)
  }
}


calculateTopicity<-function(file.list, text.names, n.div, num.sample=NA){
  library(dplyr)
  all.texts<-lapply(file.list, function(x) scan(x, what='character', sep="\n"))
  all.texts<-lapply(all.texts, function(x) paste(x, what=" "))
  all.texts<-lapply(all.texts, function(x) cleanText(x))
  all.texts<-mapply(function(x,y) divideTexts(x, n.div, y), all.texts, text.names, SIMPLIFY=F)
  if(!is.na(num.sample)){
    all.texts<-lapply(all.texts, function(x) sampleText(x, num.sample)) 
  }
  bad.index<-which(is.na(all.texts))
  if(length(bad.index)>0){
    all.texts<-all.texts[-bad.index]
    text.names<-text.names[-bad.index]
  }
  topic.sheets<-mapply(function(x,y) topicity(x,y), all.texts, text.names, SIMPLIFY=F)
  full.topicity<-bind_rows(topic.sheets)
  return(full.topicity)
}