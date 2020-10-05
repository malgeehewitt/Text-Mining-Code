cleanText<-function(raw.text, dictionary=NULL){
  split.text<-unlist(strsplit(raw.text, ""))
  split.text<-tolower(split.text)
  split.text<-split.text[which(split.text %in% c(letters, LETTERS, " "))]
  split.text<-paste(split.text, collapse="")
  split.text<-unlist(strsplit(split.text, " "))
  if(!is.null(dictionary)){
    split.text<-split.text[which(split.text %in% dictionary)]
  }
  return(split.text)
}

textWordCount<-function(file.name, dictionary){
  raw.text<-scan(file.name, what='character', sep="\n", quiet=T)
  raw.text<-paste(raw.text, collapse=" ")
  clean.text<-cleanText(raw.text, dictionary)
  clean.table<-table(clean.text)
  return(clean.table)
}

corpusWordCount<-function(folder, dictionary){
  all.files<-list.files(folder, full.names=T)
  all.word.tables<-lapply(all.files, function(x) textWordCount(x, dictionary=dictionary))
  all.word.tables<-unlist(all.word.tables)
  master.table<-tapply(all.word.tables, names(all.word.tables), sum)
  return(master.table)
}


gatherPairsText<-function(text.filename, dictionary){
  raw.text<-scan(text.filename, what='character', sep="\n", quiet=T)
  raw.text<-paste(raw.text, collapse=" ")
  clean.text<-cleanText(raw.text)
  all.counts<-table(clean.text)
  first.index<-seq(1, length(clean.text)-1, by=1)
  second.index<-first.index+1
  first.index<-clean.text[first.index]
  second.index<-clean.text[second.index]
  bigram.table<-data.frame(first.index, second.index, stringsAsFactors = F)
  if(!is.null(dictionary)){
    all.counts<-all.counts[which(names(all.counts) %in% dictionary)]
    bigram.table<-bigram.table[which(bigram.table[,1] %in% dictionary),]
    bigram.table<-bigram.table[which(bigram.table[,2] %in% dictionary),]
  }
  bigram.pairs<-paste(bigram.table[,1], bigram.table[,2], sep="_")
  anti.bigram.pairs<-paste(bigram.table[,2], bigram.table[,1], sep="_")
  all.bigrams<-c(bigram.pairs, anti.bigram.pairs)
  all.bigrams<-table(all.bigrams)
  return(list(all.bigrams, all.counts))
}


calculatePMI<-function(bigram, bigram.name, total.bigrams, word.counts){
  scaled.bigram<-bigram/total.bigrams
  bigram.names<-unlist(strsplit(bigram.name, "_"))
  bigram.counts<-word.counts[which(names(word.counts) %in% bigram.names)]
  pmi.score<-log(scaled.bigram/prod(bigram.counts))
  return(pmi.score)
}

createSVDVector<-function(pmi.dist.table, column.cut=10){
  svd.model<-svd(pmi.dist.table)
  svd.model<-svd.model$u
  svd.model<-svd.model[,1:column.cut]
  rownames(svd.model)<-rownames(pmi.dist.table)
  return(svd.model)
}

populateMatrix <- function(bigram.vector){
  unique.names <- unlist(lapply(names(bigram.vector), function(x) unlist(strsplit(x, "_"))))
  unique.names <- unique(unique.names)
  unique.names <- sort(unique.names)
  pmi.matrix <- matrix(rep(0, (length(unique.names)*length(unique.names))), nrow = length(unique.names))
  for (i in 1:length(bigram.vector)){
    pmi <- bigram.vector[i]
    cur.name <- names(bigram.vector[i])
    name.split <- unlist(strsplit(cur.name, "_"))
    pmi.matrix[which(rownames(pmi.matrix) == name.split[1], which(colnames(pmi.matrix) == name.split[2]))] <- pmi
  }
  return (pmi.matrix)
}

populateMatrixApply<-function(bigram.pmis){
  bigram.names<-names(bigram.pmis)
  unique.names<-unique(unlist(lapply(bigram.names, function(x) unlist(strsplit(x, "_")))))
  first.element<-rep(unique.names, length(unique.names))
  second.element<-sort(first.element)
  possible.bigrams<-paste(first.element, second.element, sep="_")
  remove(first.element)
  remove(second.element)
  all.pmi<-rep(1, length(possible.bigrams))
  names(all.pmi)<-possible.bigrams
  remove(possible.bigrams)
  all.pmi<-c(all.pmi, bigram.pmis)
  all.pmi<-tapply(all.pmi, names(all.pmi), sum)
  print(length(all.pmi))
  pmi.matrix<-matrix(all.pmi, ncol=length(unique.names), byrow=F)
  rownames(pmi.matrix)<-unique.names
  colnames(pmi.matrix)<-unique.names
 return(pmi.matrix)
}

#primary function
#corpus.folder is the folder of texts that you want your code to parse
#output is the name of the FILE that you want to output with the model
#dictionary is used to cut the list of words down to a manageable size by only retaining words in a dictionary file

pmiVector<-function(corpus.folder, output, dictionary){
  ptm<-proc.time()
  all.filenames<-list.files(corpus.folder, full.names=T)
  library(parallel)
  #n.core<-detectCores()-4
  #cluster.proc<-makeCluster(n.core, type="FORK")
  #all.bigram.lists<-parLapply(cluster.proc, all.filenames, function(x) gatherPairsText(x, dictionary))
  #stopCluster(cluster.proc)
  all.bigram.lists<-lapply(all.filenames, function(x) gatherPairsText(x, dictionary))
  all.bigrams<-unlist(lapply(all.bigram.lists, function(x) x[1]))
  all.word.counts<-unlist(lapply(all.bigram.lists, function(x) x[2]))
  remove(all.bigram.lists)
  all.bigrams<-tapply(all.bigrams, names(all.bigrams), sum)
  all.word.counts<-tapply(all.word.counts, names(all.word.counts), sum)
  all.bigrams<-all.bigrams[order(names(all.bigrams))]
  all.word.counts<-all.word.counts[order(names(all.word.counts))]
  total.bigrams<-sum(all.word.counts)-1
  all.word.counts<-all.word.counts/sum(all.word.counts)
  bigram.names <- names(all.bigrams)
  bigram.pmis<-unlist(mapply(function(x,y) calculatePMI(x, y, total.bigrams, all.word.counts), all.bigrams, bigram.names, SIMPLIFY = F))
  names(bigram.pmis)<-names(all.bigrams)
  #return (bigram.pmis)
  remove(all.bigrams)
  remove(all.word.counts)
  #bigram.names<-unlist(lapply(names(bigram.pmis), function(x) unlist(strsplit(x, "_"))))
  #unique.words<-unique(bigram.names)
  #bigram.matrix<-matrix(bigram.pmis, ncol=length(unique.words), byrow=T)
  #rownames(bigram.matrix)<-unique.words
  #colnames(bigram.matrix<-unique.words)
  print(proc.time()-ptm)
  ptm<-proc.time()
  bigram.matrix <- populateMatrixApply(bigram.pmis)
  print(proc.time()-ptm)
  ptm<-proc.time()
  bigram.svd<-createSVDVector(bigram.matrix, column.cut=20)
  print(proc.time()-ptm)
  write.csv(bigram.svd, file=output)
  return(bigram.svd)
}

  
  
