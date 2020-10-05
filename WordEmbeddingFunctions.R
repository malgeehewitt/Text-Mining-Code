#the code in this file
#1) allows for training a vector model on a folder of text files. Depending on the desired model, gloveModelVec will train an embedding model using GloVe and w2vModelVec will train
  #an embedding model using Mikolov's word2vec
#2)Cleans the model, taking out stopwords, numbers and non-dictionary words
#3)For any model, allows users to gather the top associated terms with a word (getTopTerms) or do math along vectors (vectorMath)
#4)Vizualize the vectors as a tSNE plot (vizualizeVectors) or as a spectrum along a single vector (vectorPlot)

#All plotting and cleaning should work with any vector model

#Requires the packages:
#text2vec
#wordVectors
#magrittr
#rTnse

#TO TRAIN A Model:
#Glove: my.model<-gloveModelVec(my.data.source, g.dimension=number of dimensions)
#word2vec: my.model<-w2vModelVec(my.data.source, g.dimension=number of dimensions)

#TO CLEAN A Model:
#dict<-scan("dictionaryfilename.txt", what='character', sep="\n)
#my.model<-cleanModel(my.model, dictionary=dict)

#To get terms
#getTopTerms(my.model, "word I want")

#To do Vector Math
#vectorMath("equation", my.model)



#GloVe word vector from texts: input (text.source) is indicated by source.type - a directory ("dir") or a vector of texts ("vector")
#stop.words removes standard stopwords from the tm package
#skip.gram indicates an integer ("L") for the skip.gram size
#term_min indicates the minimum number of times a word has to be in the model to be included
#g.dimension is the number of dimensions the model will have at the end (used for comparing vectors)
gloveModelVec<-function(text.source, stop.words=F, skip.gram=5L, term_min=5L, g.dimension=50, source.type="dir"){
  library(text2vec)
  library(magrittr)
  if(source.type=="dir"){
    texts<-lapply(list.files(text.source, full.names=T), function(x) scan(x, what='character', sep="\n", quiet=T))
    texts<-unlist(lapply(texts, function(x) paste(x, collapse=" ")))
    all.corp<-texts
  } else {
    all.corp<-text.source
  }
  tokens<-all.corp %>% tolower %>% word_tokenizer
  it<-itoken(tokens, progressbar = F)
  if(stop.words){
    library(tm)
    stop.word.list<-stopwords("en")
    vocab<-create_vocabulary(it, stopwords=stop.word.list)
  } else {
    vocab<-create_vocabulary(it)
  }
  vocab<-create_vocabulary(it)
  vocab<-prune_vocabulary(vocab, term_count_min=term_min)
  vectorizer<-vocab_vectorizer(vocab)
  tcm<-create_tcm(it, vectorizer, skip_grams_window=skip.gram)
  #glove.model<-GloVe$new(word_vectors_size=g.dimension, vocabulary=vocab, x_max=10)
  glove.model<-GloVe$new(g.dimension, x_max=10)
  word.vectors.main<-glove.model$fit_transform(tcm, n_iter=20)
  word.vectors.context<-glove.model$components
  word.vectors<-word.vectors.main + t(word.vectors.context)
  return(word.vectors)
}

#similar to above but implements Ben Schmidt's R package for Mikolov's word2vec in C
w2vModelVec<-function(text.source, skip.gram, iter, g.dimension, temp.output, model.output){
  library(wordVectors)
  library(magrittr)
  prep_word2vec(origin=text.source, destination=temp.output, lowercase=T, n.threads=1)
  word2vec.model<-train_word2vec(temp.output, model.output, vectors=g.dimension, threads=n.threads, iter=iter, window=skip.gram)
  return(word2vec.model)
}

#depreciated function for paralellizing vectors - requires a hash that is not yet implemented
gloveModelIter<-function(source.dir, min.terms= 5L, skip.gram=5L, g.dimension=50, num.iter=20, nworkers=4){
  library(text2vec)
  library(doParallel)
  library(magrittr)
  registerDoParallel(nworkers)
  all.files<-list.files(source.dir, full.names = T)
  print(length(all.files))
  jobs<-all.files %>% split_into(nworkers) %>% lapply(function(x) x %>% ifiles()) %>% lapply(itoken, tolower, word_tokenizer)
  #itoken(progressbar=F, preprocessor=tolower, tokenizer=word_tokenizer) 
  #dir.iter<-idir(source.dir)
  #it<-itoken(dir.iter, preprocessor = tolower, tokenizer=word_tokenizer)
  #vocab<-create_vocabulary(jobs)
  vocab<-create_vocabulary(jobs)
  vocab<-prune_vocabulary(vocab, term_count_min = min.terms)
  #vectorizer<-vocab_vectorizer(vocab, grow_dtm=F, skip_grams_window=skip.gram)
  vectorizer<-vocab_vectorizer(vocab)
  tcm<-create_tcm(jobs, vectorizer)
  glove <- GlobalVectors$new(word_vectors_size=g.dimension, vocabulary = vocab, x_max=10)
  glove$fit(tcm, n_iter=num.iter)
  word_vectors<-glove$get_word_vectors()
  #write.csv(word_vectors, file="EccoWordVectors.csv")
  return(word_vectors)
}

#remove paratext from a complete vector model
#glove.model is the output of gloveModelVec
#numbers indicates whether to remove words that resolve to numberic types (1,2,3,4,etc)
#roman characters indicates whether to remove the first twenty lower-case roman numerals (typically page numbers)
#letters indicates whether to remove single letters
#dictionary offers the opportunity to filter the model thorugh a vector of dictionary words (removes ocr errors and proper nouns)
cleanModel<-function(glove.model, numbers=T, roman.characters=T, letters=T, dictionary=NULL){
  terms<-rownames(glove.model)
  scrub.index<-NULL
  if(numbers){
    scrub.index<-c(scrub.index, suppressWarnings(which(!is.na(as.numeric(terms)))))
  }
  if(roman.characters){
    scrub.index<-c(scrub.index, which(terms %in% c('i','ii','iii', 'iv', 'v', 'vi', 'vii', 'viii', 'ix', 'x', 'xi', 'xiii', 'xiv', 'xv', 'xvi', 'xvii', 'xviii', 'xix', 'xx')))
  }
  if(letters){
    scrub.index<-c(scrub.index, which(terms %in% c(letters, LETTERS)))
  }
  if(!is.null(scrub.index)){
    glove.model<-glove.model[-scrub.index,]
  }
  if(!is.null(dictionary)){
    glove.model<-glove.model[which(rownames(glove.model) %in% dictionary),]
  }
  return(glove.model)
}

#code to import a pre-trained vector model from a file
importVectorModel<-function(model.filename){
  vector.model<-read.csv(model.filename, header=T, stringsAsFactors=F, row.names=1)
  vector.model<-as.matrix(vector.model)
  return(vector.model)
}


#retrieve the top n.return words from a model closest to term - can be used as a filtering wordlist in plotting functions
#gloveModel is the output of gloveModelVec

getTopTerms<-function(gloveModel, term, n.return=25, method="cosine"){
  library(text2vec)
  if(length(term)==1){
    term.results<-gloveModel[term, ,drop=F]
  } else {
    term.results<-term
  }
  if(method=="cosine"){
    vector.results<-sim2(x=gloveModel, y=term.results, method='cosine', norm='l2')
    vector.results<-head(sort(vector.results[,1], decreasing=T),n.return)
  } else {
    vector.results<-apply(gloveModel, 1, function(x) dist(rbind(term.results, x), method=method))
    vector.results<-head(sort(vector.results), n.return)
  }
  return(vector.results)
}

#given term1 and term2, and a glove model, subtracts the vector for term1 from the vector for term2 
#returns the top n.return items from the resulting vector
subtractVec<-function(term1, term2, gloveModel, n.return=25){
  if(length(term1)==1){
    term1<-gloveModel[term1, , drop=F]
  }
  if(length(term2)==1){
    term2<-gloveModel[term2, , drop=F]
  }
  result.vector<-term1-term2
  result.vector<-getTopTerms(gloveModel, result.vector, n.return)
  return(result.vector)
}

#given term1 and term2, and a glove model, adds the vector for term1 to the vector for term2 
#returns the top n.return items from the resulting vector
addVec<-function(term1, term2, gloveModel, n.return=25){
  if(length(term1)==1){
    term1<-gloveModel[term1, , drop=F]
  }
  if(length(term2)==1){
    term2<-gloveModel[term2, , drop=F]
  }
  result.vector<-term1+term2
  result.vector<-getTopTerms(gloveModel, result.vector, n.return)
  return(result.vector)
}

#given term1 and term2, and a glove model, multiplies the vector for term1 by the vector for term2 
#returns the top n.return items from the resulting vector
multiplyVec<-function(term1, term2, gloveModel, n.return=25){
  if(length(term1)==1){
    term1<-gloveModel[term1, , drop=F]
  }
  if(length(term2)==1){
    term2<-gloveModel[term2, , drop=F]
  }
  result.vector<-term1*term2
  result.vector<-getTopTerms(gloveModel, result.vector, n.return)
  return(result.vector)
}

#given term1 and term2, and a glove model, divides the vector for term1 by the vector for term2 
#returns the top n.return items from the resulting vector
divideVec<-function(term1, term2, gloveModel, n.return=25){
  if(length(term1)==1){
    term1<-gloveModel[term1, , drop=F]
  }
  if(length(term2)==1){
    term2<-gloveModel[term2, , drop=F]
  }
  result.vector<-term1/term2
  result.vector<-getTopTerms(gloveModel, result.vector, n.return)
  return(result.vector)
}

#wrapper for compund math operations
#vector formula is in the form of term operation term operation (e.g. term1+term2-term3)
#glove.model is the result of gloveModelVec
#function returns top n.return terms from the resulting vector
#NOTE: the function words in order from left to right rather than following the order of operations
vectorMath<-function(vector.formula, glove.model, n.return=25, method="cosine"){
  operation.signs<-c("*", "/", "+", "-")
  chars<-unlist(strsplit(vector.formula, ''))
  operation.index<-which(chars %in% operation.signs)
  operations<-chars[operation.index]
  start.words<-c(1, (operation.index+1))
  end.words<-c((operation.index-1), length(chars))
  words<-unlist(mapply(function(x,y) paste(chars[x:y], collapse=''), start.words, end.words))
  word.vectors<-lapply(words, function(x) glove.model[x, ,drop=F])
  composite.formula<-as.list(rep(NA, length(word.vectors)+length(operations)))
  composite.formula[seq(1, length(composite.formula), by=2)]<-word.vectors
  composite.formula[seq(2, length(composite.formula), by=2)]<-operations
  sign.index<-which(composite.formula %in% operation.signs)
  for(i in 1:length(sign.index)){
    curr.index<-sign.index[i]
    op.sign<-composite.formula[[curr.index]]
    vec1<-composite.formula[[(curr.index-1)]]
    vec2<-composite.formula[[(curr.index+1)]]
    if(op.sign=="*"){
      new.vec<-vec1*vec2
    } else if(op.sign=="/"){
      new.vect<-vec1/vec2
    } else if(op.sign=="+"){
      new.vec<-vec1+vec2
    } else if(op.sign=="-"){
      new.vec<-vec1-vec2
    }
    composite.formula[(curr.index+1)]<-list(new.vec)
    composite.formula<-composite.formula[-c((curr.index-1), curr.index)]
    sign.index<-sign.index-2
  }
  result.vector<-composite.formula[[1]]
  result.vector<-getTopTerms(glove.model, result.vector, n.return, method=method)
  return(result.vector)
}

#create a pdf file of a plot of all words in a glove model, colored for kmeans groups
#glove.model is the result of gloveModelVec
#filename is the output filename for the .pdf file
#k.groups indicates how many k-means groups are modeled
#pdf.height is the height in inches of the resulting .pdf file
#pdf.width is the width in inches of the resulting .pdf file
#wordlist filters the final plot for only words in the wordlist (used in conjunction with getTopTerms it can filter results)
#NOTE: Rtsne is a slow but accurate dimension reduction - on vectors with more than 10000 terms, it can take upwards of an hour
visualizeVectors<-function(glove.model, filename, k.groups, pdf.height=50, pdf.width=75, wordlist=NULL, kCenters=F){
  library(Rtsne)
  if(!is.null(wordlist)){
    glove.model<-glove.model[which(rownames(glove.model) %in% wordlist),]
  }
  if(ncol(glove.model)>2){
    reduced.model<-Rtsne(glove.model, check_duplicates=F)
    reduced.matrix<-reduced.model$Y
    rownames(reduced.matrix)<-rownames(glove.model)
  } else {
    reduced.matrix<-glove.model
  }
  terms<-rownames(glove.model)
  k.means<-kmeans(reduced.matrix, k.groups)
  k.centers<-k.means$centers
  k.clusters<-factor(k.means$cluster)
  plot.table<-data.frame(reduced.matrix, terms, k.clusters)
  colnames(plot.table)<-c("Tsne1", "Tsne2", "Term", "Cluster")
  if(kCenters){
    k.center.addition<-data.frame(k.centers, rep("KCenter", nrow(k.centers)), rep(as.character(k.groups+1), nrow(k.centers)))
    colnames(k.center.addition)<-c("Tsne1", "Tsne2", "Term", "Cluster")
    plot.table<-rbind(plot.table, k.center.addition)
  }
  colnames(plot.table)<-c("Tsne1", "Tsne2", "Term", "Cluster")
  library(ggplot2)
  embed.plot<-ggplot(plot.table, aes(x=Tsne1, y=Tsne2, color=Cluster, label=Term))+geom_text(size=2)#+geom_point()+geom_text(size=0.7, color="black")
  pdf(filename, height=pdf.height, width = pdf.width)
  print(embed.plot)
  dev.off()
  return(plot.table)
}

#plot terms in a fector model that lie on an axis between term1 and term2 (both words in the model)
#this reoirients the plot of terms so that, from left to right, words become less associated with term1 and more with term 2
#vec.model is the result of gloveModelVec
#filename is the output filename for the .pdf file
#pdf.height is the height in inches of the resulting .pdf file
#pdf.width is the width in inches of the resulting .pdf file
#wordlist filters the final plot for only words in the wordlist (used in conjunction with getTopTerms it can filter results)
vectorPlot<-function(term1, term2, vec.model, filename, pdf.height=50, pdf.width=75, wordlist=NULL){
  library(ggplot2)
  #vec.model<-vec.model[which(rownames(vec.model) %in% wordlist),]
  subtract.results<-subtractVec(term1, term2, vec.model, nrow(vec.model))
  wordlist<-c(wordlist, term1, term2)
  subtract.results<-subtract.results[which(names(subtract.results) %in% wordlist)]
  terms<-names(subtract.results)
  x.values<-as.vector(subtract.results)
  y.values<-sample(seq(-1,1,by=0.0001), length(x.values))
  plot.table<-data.frame(x.values, y.values, terms)
  term.plot<-ggplot(plot.table, aes(x=x.values, y=y.values, label=terms))+geom_text(color="black", size=4)+ggtitle(paste("Vector of:", term1, "to", term2, sep=" "))
  pdf(filename, height=pdf.height, width=pdf.width)
  print(term.plot)
  dev.off()
  return(plot.table)
}