#function for finding distinctive words from a corpus of texts with group assignments
#function takes:
#corpus: either as a variable of class corpus from TM package or feature table of texts by words with raw obs values
#groups: a vector of group assignments of same length or nrow as corpus
#alpha: alpha value against which to test significance
#class: a boolean flag indicating if the results are to be used for classification (if T, will limit the number of variables generated to prevent overfitting of model)
#greater.only: a boolean flag indicating if all significant values (more or less than expected) are returned or only variables appearing significantly more often
#scaling.vector: option to include a separate scaling vector- if NULL defaults to the number of words in each text
#exclude.zero: deterimines in single text or group words are retained; can be "None", "Group" or "Text" - if "None", returns all words, if "Text", excludes all words that appear in only one text; if "Group" excludes words appearing in only one group
#function returns a maxtrix of significant results, including word list, raw number of observations, observed/expected for each word, pvalue and group associated with word

getGroupMDWs<-function(corpus, groups, alpha=0.05, class=F, greater.only=T, scaling.vector=NULL, exclude.zero="None", dictionary=NULL){
  library(tm)
  
  #determine if corpus variable is a raw corpus or a feature table
  if(is.null(nrow(corpus))){
    dtm.flag<-FALSE
  } else {
    dtm.flag<-TRUE
  }
  #convert grouping variable to factor and find unique groups
  groups<-as.factor(groups)
  group.names<-levels(groups)
  
  #clean and create feature table if corpus variable is raw corpus
  if (!dtm.flag){
    print("Cleaning Texts")
    corpus<-prepCorpus(corpus)
    corpus<-hardClean(corpus)#, dictionary)
    obs.dtm<-makeDtm(corpus, stop.words=F)
    print(dim(obs.dtm))
    
  } else {
    obs.dtm<-corpus
  }

  #exclude words that only appear in one text
  if (exclude.zero=="Text"){
    num.zeros<-apply(obs.dtm, 2, function(x) length(x[x==0]))
    single.text.index<-which(num.zeros==((nrow(obs.dtm))-1))
    obs.dtm<-obs.dtm[,-single.text.index]
  }
  
  #scale dtm to find expected values per text for full corpus
  word.totals.words<-colSums(obs.dtm)
  if (is.null(scaling.vector)){
    if(!dtm.flag){
      full.dtm<-makeDtm(corpus)
    } else {
      full.dtm<-corpus
    }
    total.words<-sum(full.dtm)
    scaling.vector<-rowSums(full.dtm)
  } else {
    total.words<-sum(scaling.vector)
  }
  word.freq.corpus<-word.totals.words/total.words
  
  #For classification - determine maximum number of variables
  if (class){
    num.groups<-length(group.names)
    max.var<-nrow(obs.dtm) %/% 3
    max.var.group<-max.var %/% num.groups
  }
  
  #intitalize MDW table
  MDW.table<-NULL
  
  #for each unique group, extract corresponding texts from feature table for observed values
  #and calculate expected values based on size (number of words) of group
  for (i in 1:length(group.names)){
    print(group.names[i])
    group.index<-which(groups == group.names[i])
    group.texts<-obs.dtm[group.index,]
    if(length(group.index)==1){
      group.texts<-t(as.matrix(group.texts))
    }
    obs.values<-colSums(group.texts)
    group.text.lengths<-rowSums(group.texts)
    total.group.words<-sum(scaling.vector[group.index])
    exp.values<-word.freq.corpus*total.group.words
    names(obs.values)<-colnames(obs.dtm)
    word.probs<-word.freq.corpus

    #exclude words that only appear in one group
    if (exclude.zero=="Group"){
      non.group.freq<-colSums(obs.dtm[-group.index,])
      zero.index<-which(non.group.freq==0)
      if (length(zero.index>0)){
        obs.values<-obs.values[-zero.index]
        word.probs<-word.probs[-zero.index]
        exp.values<-exp.values[-zero.index]
      }
    }
    #end exclusion
    
    #greater.only flag - only retain words with frequencies greater than expected
    if (greater.only){
      lessThan.index<-which(exp.values>obs.values)
      obs.values<-obs.values[-lessThan.index]
      word.probs<-word.probs[-lessThan.index]
    }
    
    #get mdws for single unique group variable
    group.MDWs<-getMDWs(obs.values, word.probs, total.group.words, alpha)
    if(!is.null(group.MDWs)){
      #if only one MDW is returned, create matrix out of vector
      if(is.null(ncol(group.MDWs))){
        group.MDWs<-as.matrix(group.MDWs)
      }
      #associate group name with group MDWs
      MDW.group<-rep(group.names[i], ncol(group.MDWs))
      group.MDWs<-rbind(group.MDWs, MDW.group)
      #for classification - reduce number of variables to avoid overfitting
      if (class){
        group.MDWs<-group.MDWs[,1:max.var.group]
      }
      #add current group MDWs to overall MDW table
      if (i==1){
        MDW.table<-group.MDWs
      } else {
        if(nrow(group.MDWs) != nrow(MDW.table)){
          print(group.MDWs)
        }
        MDW.table<-cbind(MDW.table, group.MDWs)
      }
    }
    print(ncol(group.MDWs))
  }
  
  rownames(MDW.table)<-c("Word", "N_Obs", "Obs_Exp", "P Value", "Sig_Group")
  colnames(MDW.table)<-NULL
  MDW.table<-MDW.table[,order(-as.numeric(MDW.table[3,]))]
  MDWs<-MDW.table[1,]
  #resolve duplicated MDWs by retaining the group incidence with the highest Obs/Exp
  duplicate.index<-which(duplicated(MDWs))
  if (length(duplicate.index)>0){
    MDW.table<-MDW.table[,-duplicate.index]
  }
  
  return(t(MDW.table))
}

#tm package function for corpus preparation
prepCorpus<-function(corpus){
  library(tm)
  corpus<-tm_map(corpus, content_transformer(removePunctuation))
  corpus<-tm_map(corpus, content_transformer(removeNumbers))
  return(corpus)
}

#text clean - remove all non-letter, space and EoL characters
hardClean<-function(corpus){#, dict.list){
  corp.names<-names(corpus)
  text.vector<-NULL
  library(tm)
  corpus<-tm_map(corpus, content_transformer(tolower))
  for (i in 1:length(corpus)){
    #curr.text<-corpus[[i]]$content
    curr.text<-corpus[[i]]$content
    curr.text<-unlist(strsplit(curr.text, ' '))
    curr.text<-paste(curr.text, collapse=' ')
    curr.text.split<-unlist(strsplit(curr.text, ''))
    dash.index<-which(curr.text.split == "-")
    curr.text.split[dash.index]<-' '
    bad.index<-which(!(curr.text.split %in% letters))
    space.index<-which(bad.index %in% which(curr.text.split == ' '))
    if (length(space.index) < length(bad.index)){
      bad.index<-bad.index[-space.index]
      curr.text.clean<-curr.text.split[-bad.index]
    } else {
      curr.text.clean<-curr.text.split
    }
    curr.text.clean<-paste(curr.text.clean, collapse='')
    #keep only dictionary words
    # curr.text.clean<-unlist(strsplit(curr.text.clean, " "))
    # curr.text.clean<-curr.text.clean[which(curr.text.clean %in% dict.list)]
    # curr.text.clean<-paste(curr.text.clean, collapse=" ")
    #keep only dictionary words
    text.vector<-c(text.vector, curr.text.clean)
  }
  text.vector.source<-VectorSource(text.vector)
  new.corpus<-Corpus(text.vector.source)
  #names(new.corpus)<-corp.names
  return (new.corpus)
}

#tm function-create a scaled DTM as a matrix
makeScaledDtm<-function(corpus, stop.words=F, sparsity=1){
  library(tm)
  if (stop.words == T){
    params<-list(stopwords=TRUE, wordLengths=c(1,Inf))
  } else {
    params<-list(wordLengths=c(1,Inf))
  }
  full.dtm<-DocumentTermMatrix(corpus)
  full.dtm<-as.matrix(full.dtm)
  scaling<-rowSums(full.dtm)
  corpus.dtm<-DocumentTermMatrix(corpus, control=params)
  if (sparsity!=1){
    corpus.dtm<-removeSparseTerms(corpus.dtm, sparsity)
  }
  corpus.dtm<-as.matrix(corpus.dtm)
  dtm.scaled<-corpus.dtm[,1:ncol(corpus.dtm)]/scaling
  return(dtm.scaled)
}

#tm function - make unscaled DTM (using raw Obs values)
makeDtm<-function(corpus, stop.words=F, sparsity=.7){
  library(tm)
  if (stop.words == T){
    params<-list(stopwords=TRUE, wordLengths=c(1,Inf))
  } else {
    params<-list(wordLengths=c(1,Inf))
  }
  corpus.dtm<-DocumentTermMatrix(corpus, control=params)
  if (sparsity!=1){
    corpus.dtm<-removeSparseTerms(corpus.dtm, sparsity)
  }
  corpus.dtm<-as.matrix(corpus.dtm)
  return(corpus.dtm)
}

#function for creating grouping variables - not run as part of primary function
makeGroups<-function(corpus, group.part=1, num.parts=2){
  text.names<-names(corpus)
  text.names<-unlist(strsplit(text.names, ".txt"))
  group.names<-unlist(strsplit(text.names, "_"))
  group.index<-seq(group.part,length(group.names), by=num.parts)
  group.names<-group.names[group.index]
  return(group.names)
}


#given a vector of obs values for a single group, a vector of frequencies for the same group,
#a total size for the group and an alpha value, creates a 2x2 contingency table of presence/absence
#per word and uses a fishers exact text to compare the p value of each contingency table to alpha
#returns matrix of values per sigificant word (word, nobs, obs/exp and pvalue) or NULL if
#no signifcant words are found
getMDWs<-function(obs, freq, total.group.size, alpha){
  MDW.table<-NULL
  words<-names(obs)
  for (j in 1:length(obs)){
    #exp<-freq[j]*total.group.size
    #exp<-as.integer(exp)
    present<-c(obs[j], as.integer(freq[j]*total.group.size))
    absent<-total.group.size-present
    contingency.table<-rbind(present,absent)
    #print(contingency.table)
    word.p<-fisher.test(contingency.table, alternative="g")$p.value
    if (word.p<alpha){
      curr.MDW<-c(words[j], obs[j], obs[j]/(freq[j]*total.group.size), word.p)
      MDW.table<-cbind(MDW.table, curr.MDW)
    } 
  }
  MDW.table.sort<-MDW.table[,order(-as.numeric(MDW.table[3,]))]
  return(MDW.table.sort)
}