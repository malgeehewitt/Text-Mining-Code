text.split<-function(curr.text, num.div){
  text.divs<-NULL
  curr.length<-length(curr.text)
  div.length<-curr.length%/%num.div
  for (i in 1:num.div){
    curr.div<-curr.text[(((i-1)*div.length)+1):(i*div.length)]
    text.divs<-rbind(text.divs, curr.div)
  }
  return (text.divs)
}

field.count<-function(chopped.text, field){
  div.counts<-NULL
  for (i in 1:nrow(chopped.text)){
    n.hits<-length(which(chopped.text[i,] %in% field))
    div.counts<-c(div.counts, n.hits)
  }
  return(div.counts)
}

count.field<-function(text.divisions, field){
  hit.count<-NULL
  for (i in 1:nrow(text.divisions)){
    num.hits<-length(which(text.divisions[i,] %in% field))
    num.hits<-num.hits/ncol(text.divisions)
    curr.hits<-c(i, num.hits)
    hit.count<-rbind(hit.count, curr.hits)
  }
  colnames(hit.count)<-c("section", "percent")
  return (hit.count)
}
 
#given a single text, a number of divisions and a single list of words returns a ggplot of the percentage of words from the
#list within each division  
plot.fields<-function(corpus.text, division.length, field.words, poly.degree=8){
  corpus.words<-strsplit(corpus.text, " ")
  corpus.words<-unlist(corpus.words)
  corpus.parts<-text.split(corpus.words, division.length)
  corpus.fields<-count.field(corpus.parts, field.words)
  corpus.fields<-as.data.frame(corpus.fields)
  p<-ggplot(corpus.fields, aes(x=section, y=percent, color=percent)) + geom_point()+geom_line()
  #p1<-p + stat_smooth(method="lm", formula=y~poly(x, poly.degree), se=F)
  return (corpus.fields)
}


#takes a corpus text, divides it into n.division parts and returns a new corpus of its parts
textChop<-function(corpus.text, n.division){
  text.vector<-NULL
  library(tm)
  corpus.words<-strsplit(corpus.text, " ")
  corpus.words<-unlist(corpus.words)
  corpus.parts<-text.split(corpus.words, n.division)
  for (i in 1:nrow(corpus.parts)){
    text.part<-paste(corpus.parts[i,], collapse=' ')
    text.vector=c(text.vector, text.part)
  }
  text.vectorSource<-VectorSource(text.vector)
  chopped.corpus<-Corpus(text.vectorSource, readerControl=list(language="English"))
  return(chopped.corpus)
}

#given a corpus and a part length divides each text in the corpus into parts of division.size length
#and returns all parts as a new corpus
textChopStatic<-function(corpus, division.size){
  text.vector<-NULL
  split.text.names<-NULL
  names.texts<-names(corpus)
  library(tm)
  num_text<-length(corpus)
  ver<-R.Version()
  ver.num<-as.numeric(ver$minor)
  length.table<-NULL
  start.indicies<-NULL
  end.indicies<-NULL
  for (i in 1:num_text){
    #print(names.texts[i])
    if(ver.num>=1.1){
      curr.text<-corpus[[i]]$content
    } else {
      curr.text<-corpus[[i]]
    }
    curr.name=names.texts[i]
    text.words<-unlist(strsplit(curr.text, " "))
    text.size<-length(text.words)
    num.parts<-text.size %/% division.size
    #print(num.parts)
    start.index<-0
    #print(num.parts)
    for (j in 0:num.parts){
      if (j == num.parts){
        if(start.index<length(text.words)){
          curr.chunk=text.words[start.index:length(text.words)]
        }
      } else {
        end.index<-start.index+division.size
        curr.chunk=text.words[start.index:end.index]
      }
      start.index<-end.index+1
      curr.chunk.name<-c(curr.name, j)
      curr.chunk.name<-paste(curr.chunk.name, collapse="_")
      curr.chunk=paste(curr.chunk, collapse=" ")
      text.vector<-c(text.vector, curr.chunk)
      split.text.names<-c(split.text.names, curr.chunk.name)
    }
  }
  split.vector<-VectorSource(text.vector)
  split.corpus<-Corpus(split.vector, readerControl=list(language="English"))
  names(split.corpus)<-split.text.names
  return(split.corpus)
}
  
  

plot.field<-function(corpus.text, division.length, field.words){
  corpus.words<-strsplit(corpus.text, " ")
  corpus.words<-unlist(corpus.words)
  corpus.parts<-text.split(corpus.words, division.length)
  corpus.fields<-count.field(corpus.parts, field.words)
  return(corpus.fields)
}

#takes a corpus text, a list of fields in a matrix and a division length and returns the data frame with the part values
multi.field.table<-function(corpus.text, division.length, fields, poly.degree=8){
  all.fields<-NULL
  aggregated<-unlist(fields)
  aggregated<-aggregated[!duplicated(aggregated)]
  sense.names<-colnames(fields)
  for (i in 1:ncol(fields)){
    curr.fields<-fields[,i]
    #bad.index<-which(curr.fields == "")
    #curr.fields<-curr.fields[-bad.index]
    curr.fields.matrix<-plot.field(corpus.text, division.length, curr.fields)
    curr.sense<-sense.names[i]
    sense.col<-rep(curr.sense, division.length)
    curr.fields.matrix<-cbind(curr.fields.matrix, sense.col)
    all.fields<-rbind(all.fields, curr.fields.matrix)
  }
  #print(aggregated)
  curr.fields.matrix<-plot.field(corpus.text, division.length, aggregated)
  curr.sense<-"aggregated"
  sense.col<-rep(curr.sense, division.length)
  curr.fields.matrix<-cbind(curr.fields.matrix, sense.col)
  all.fields<-rbind(all.fields, curr.fields.matrix)
  colnames(all.fields)<-c("division", "percent", "sense")
  division.names<-all.fields[,1]
  all.fields<-as.data.frame(all.fields, stringsAsFactor=False)
  all.fields[,3]<-as.factor(all.fields[,3])
  division.names<-as.numeric(division.names)
  all.fields[,1]<-division.names
  return(all.fields)
}

#takes a corpus text, a list of fields in a matrix and a division length and plots the fields together with their aggregated value
plot.multi.fields<-function(corpus.text, division.length, fields, poly.degree=8){
  require(ggplot2)
  all.fields<-NULL
  aggregated<-unlist(fields)
  aggregated<-aggregated[!duplicated(aggregated)]
  sense.names<-colnames(fields)
  for (i in 1:ncol(fields)){
    curr.fields<-fields[,i]
    #bad.index<-which(curr.fields == "")
    #curr.fields<-curr.fields[-bad.index]
    curr.fields.matrix<-plot.field(corpus.text, division.length, curr.fields)
    curr.sense<-sense.names[i]
    sense.col<-rep(curr.sense, division.length)
    curr.fields.matrix<-cbind(curr.fields.matrix, sense.col)
    all.fields<-rbind(all.fields, curr.fields.matrix)
  }
  #print(aggregated)
  curr.fields.matrix<-plot.field(corpus.text, division.length, aggregated)
  curr.sense<-"aggregated"
  sense.col<-rep(curr.sense, division.length)
  curr.fields.matrix<-cbind(curr.fields.matrix, sense.col)
  all.fields<-rbind(all.fields, curr.fields.matrix)
  colnames(all.fields)<-c("division", "percent", "sense")
  division.names<-all.fields[,1]
  all.fields<-as.data.frame(all.fields, stringsAsFactor=False)
  all.fields[,3]<-as.factor(all.fields[,3])
  division.names<-as.numeric(division.names)
  all.fields[,1]<-division.names
  p<-ggplot(all.fields, aes(x=division, y=percent, colour=sense, group=sense))+geom_point()#aes(shape=sense))
  p1<-p+stat_smooth(method="lm", formula=y~poly(x,poly.degree), se=F)
  return(p1)
}

#takes a corpus text, division length and multiple lists of fields in a matrix and plots them without the aggregated line
plot.multi.nonagg<-function(corpus.text, division.length, fields, poly.degree=8){
  all.fields<-NULL
  aggregated<-unlist(fields)
  aggregated<-aggregated[!duplicated(aggregated)]
  sense.names<-colnames(fields)
  for (i in 1:ncol(fields)){
    curr.fields<-fields[,i]
    #bad.index<-which(curr.fields == "")
    #curr.fields<-curr.fields[-bad.index]
    curr.fields.matrix<-plot.field(corpus.text, division.length, curr.fields)
    curr.sense<-sense.names[i]
    sense.col<-rep(curr.sense, division.length)
    curr.fields.matrix<-cbind(curr.fields.matrix, sense.col)
    all.fields<-rbind(all.fields, curr.fields.matrix)
  }
  #print(aggregated)
  #curr.fields.matrix<-plot.field(corpus.text, division.length, aggregated)
  #curr.sense<-"aggregated"
  #sense.col<-rep(curr.sense, division.length)
  #curr.fields.matrix<-cbind(curr.fields.matrix, sense.col)
  #all.fields<-rbind(all.fields, curr.fields.matrix)
  colnames(all.fields)<-c("division", "percent", "sense")
  division.names<-all.fields[,1]
  all.fields<-as.data.frame(all.fields, stringsAsFactor=False)
  all.fields[,3]<-as.factor(all.fields[,3])
  division.names<-as.numeric(division.names)
  all.fields[,1]<-division.names
  p<-ggplot(all.fields, aes(x=division, y=percent, colour=sense, group=sense))+geom_point()
  p1<-p+geom_line()#stat_smooth(method="lm", formula=y~poly(x,poly.degree), se=F)
  return(p1)
}

#takes a corpus text, division length and multiple lists of fields in a matrix and plots them without the aggregated line
#version 2 - for use with MultiPlotter.R - creates multiple plots across a number of files
external.plot.multi.nonagg<-function(corpus.text, division.length, fields, poly.degree=8, text.name){
  all.fields<-NULL
  aggregated<-unlist(fields)
  aggregated<-aggregated[!duplicated(aggregated)]
  sense.names<-colnames(fields)
  for (i in 1:ncol(fields)){
    curr.fields<-fields[,i]
    #bad.index<-which(curr.fields == "")
    #curr.fields<-curr.fields[-bad.index]
    curr.fields.matrix<-plot.field(corpus.text, division.length, curr.fields)
    curr.sense<-sense.names[i]
    sense.col<-rep(curr.sense, division.length)
    curr.fields.matrix<-cbind(curr.fields.matrix, sense.col)
    all.fields<-rbind(all.fields, curr.fields.matrix)
  }
  #print(aggregated)
  #curr.fields.matrix<-plot.field(corpus.text, division.length, aggregated)
  #curr.sense<-"aggregated"
  #sense.col<-rep(curr.sense, division.length)
  #curr.fields.matrix<-cbind(curr.fields.matrix, sense.col)
  #all.fields<-rbind(all.fields, curr.fields.matrix)
  colnames(all.fields)<-c("division", "percent", "fields")
  division.names<-all.fields[,1]
  division<-as.numeric(division.names)
  percent<-all.fields[,2]
  percent<-as.numeric(percent)
  fields<-all.fields[,3]
  fields<-factor(fields)
  all.fields<-data.frame(division, percent, fields)
  #all.fields[,3]<-as.factor(all.fields[,3])
  #division.names<-as.numeric(division.names)
  #all.fields[,1]<-division.names
  #percents<-all.fields[,2]
  #percents<-as.numeric(percents)
  #all.fields[,2]<-percents
  p<-ggplot(all.fields, aes(x=division, y=percent, colour=fields, group=fields))+geom_point(aes(shape=fields))
  p1<-p+stat_smooth(method="lm", formula=y~poly(x,poly.degree), se=F)
  p2<-p1+labs(title=text.name)+theme(axis.text=element_text(size=8))
  p3<-p2+scale_shape_manual(values=1:length(levels(all.fields$field)))
  return(p3)
}

#same as above but plots lines between points instead of a regression
plot.multi.fields.line<-function(corpus.text, division.length, fields){
  all.fields<-NULL
  aggregated<-unlist(fields)
  aggregated<-aggregated[!duplicated(aggregated)]
  sense.names<-colnames(fields)
  for (i in 1:ncol(fields)){
    curr.fields<-fields[,i]
    #bad.index<-which(curr.fields == "")
    #curr.fields<-curr.fields[-bad.index]
    curr.fields.matrix<-plot.field(corpus.text, division.length, curr.fields)
    curr.sense<-sense.names[i]
    sense.col<-rep(curr.sense, division.length)
    curr.fields.matrix<-cbind(curr.fields.matrix, sense.col)
    all.fields<-rbind(all.fields, curr.fields.matrix)
  }
  #print(aggregated)
  curr.fields.matrix<-plot.field(corpus.text, division.length, aggregated)
  curr.sense<-"aggregated"
  sense.col<-rep(curr.sense, division.length)
  curr.fields.matrix<-cbind(curr.fields.matrix, sense.col)
  all.fields<-rbind(all.fields, curr.fields.matrix)
  colnames(all.fields)<-c("division", "percent", "sense")
  division.names<-all.fields[,1]
  all.fields<-as.data.frame(all.fields, stringsAsFactor=False)
  all.fields[,3]<-as.factor(all.fields[,3])
  division.names<-as.numeric(division.names)
  all.fields[,1]<-division.names
  p<-ggplot(all.fields, aes(x=division, y=percent, colour=sense, group=sense))+geom_point()#aes(shape=sense))
  p1<-p+geom_line()
  return(p1)
}
  
 #given a corpus and a list of terms, returns a plot of the percentage of each corpus member in that list 
 plotFieldCorpus<-function(corpus.texts, fields){
  #source("GetTaggedText.r")
  hit.count<-NULL
  for (i in 1:length(corpus.texts)){
    curr.text<-corpus.texts[[i]]
    #curr.text<-striptags(curr.text)
    curr.text<-strsplit(curr.text, " ")
    curr.text<-unlist(curr.text)
    num.hits<-length(which(curr.text %in% fields))
    num.hits<-num.hits/length(curr.text)
    curr.hits<-c(i, num.hits)
    hit.count<-rbind(hit.count, curr.hits)
  }
  colnames(hit.count)<-c("section", "percent")
  corpus.fields<-hit.count
  #corpus.fields<-as.data.frame(corpus.fields)   
  #p<-ggplot(corpus.fields, aes(x=section, y=percent, color=percent)) + geom_point()
  #p1<-p + stat_smooth(method="lm", formula=y~poly(x, 8), se=F)
  return (corpus.fields)
 }
 
 #takes a corpus and multiple lists of fields in a matrix and plots them without the aggregated line
plot.corpus.multi<-function(corpus, fields, poly.degree=8){
  all.fields<-NULL
  aggregated<-unlist(fields)
  aggregated<-aggregated[!duplicated(aggregated)]
  sense.names<-colnames(fields)
  for (i in 1:ncol(fields)){
    curr.fields<-fields[,i]
    #bad.index<-which(curr.fields == "")
    #curr.fields<-curr.fields[-bad.index]
    curr.fields.matrix<-plotFieldCorpus(corpus, curr.fields)
    curr.sense<-sense.names[i]
    sense.col<-rep(curr.sense, length(corpus))
    curr.fields.matrix<-cbind(curr.fields.matrix, sense.col)
    all.fields<-rbind(all.fields, curr.fields.matrix)
  }
  #print(aggregated)
  #curr.fields.matrix<-plot.field(corpus.text, division.length, aggregated)
  #curr.sense<-"aggregated"
  #sense.col<-rep(curr.sense, division.length)
  #curr.fields.matrix<-cbind(curr.fields.matrix, sense.col)
  #all.fields<-rbind(all.fields, curr.fields.matrix)
  colnames(all.fields)<-c("segment", "percent", "field")
  division.names<-all.fields[,1]
  all.fields<-as.data.frame(all.fields, stringsAsFactor=False)
  all.fields[,3]<-as.factor(all.fields[,3])
  division.names<-as.numeric(division.names)
  all.fields[,1]<-division.names
  p<-ggplot(all.fields, aes(x=segment, y=percent, colour=field, group=field))+geom_point()#aes(shape=field))
  p1<-p+stat_smooth(method="lm", formula=y~poly(x,poly.degree), se=F)
  return(p1)
}
 
 #given a single text and a number of divisions, returns the divisions as members of a new corpus
 createSplitCorpus<-function(corpus.text, division.length){ 
  corpus.words<-strsplit(corpus.text, " ")
  corpus.words<-unlist(corpus.words)
  corpus.parts<-text.split(corpus.words, division.length)
  corpus.collapsed<-NULL
  for (i in 1:nrow(corpus.parts)){
    to.add<-paste(corpus.parts[i,], collapse=" ")
    corpus.collapsed<-rbind(corpus.collapsed, to.add)
  }
  corpus.vector<-VectorSource(corpus.collapsed)
  new.corpus<-Corpus(corpus.vector, readerControl=list(language="English"))
 }
 
plotHscore<-function(corpus.text, numdiv){
  corpustext.split<-createSplitCorpus(corpus.text, numdiv)
  corpustext.split<-tm_map(corpustext.split, removeNumbers)
  corpustext.split<-tm_map(corpustext.split, removePunctuation)
  corpustext.split.dtm<-DocumentTermMatrix(corpustext.split)
  corpustext.split.matrix<-as.matrix(corpustext.split.dtm)
  corpustext.pca<-prcomp(corpustext.split.matrix)
  plot(corpustext.pca$x[,1], corpustext.pca$x[,2])
  text(corpustext.pca$x[,1], corpustext.pca$x[,2], rownames(corpustext.split.matrix))
  return(corpustext.pca)
 } 

 plotMultiTexts<-function(corpus, numdiv, field){
  all.texts<-NULL
  text.names<-names(corpus)
  for (i in 1:length(corpus)){
    curr.text<-corpus[[i]]
    new.text.table<-plot.field(curr.text, numdiv, field)
    curr.name<-rep(text.names[i], nrow(new.text.table))
    add.table<-cbind(new.text.table, curr.name)
    all.texts<-rbind(all.texts, add.table)
  }
  colnames(all.texts)<-c("division", "percent", "texts")
  division.names<-all.texts[,1]
  all.fields<-as.data.frame(all.texts, stringsAsFactor=False)
  all.fields[,3]<-as.factor(all.fields[,3])
  division.names<-as.numeric(division.names)
  all.fields[,1]<-division.names
  p<-ggplot(all.fields, aes(x=division, y=percent, colour=texts, group=texts))+geom_point(aes(shape=texts))
  p1<-p+stat_smooth(method="lm", formula=y~poly(x,8), se=F)
  return(p1)
 }
 