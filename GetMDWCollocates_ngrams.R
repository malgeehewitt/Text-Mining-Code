collapseNgram<-function(text, ngram){
  ngram.split<-unlist(strsplit(ngram, " "))
  n<-length(ngram.split)
  text.split<-unlist(strsplit(text, " "))
  start.index<-which(text.split == ngram.split[1])
  if (length(start.index>0)){
    if (start.index[length(start.index)]+n>length(text.split)){
      start.index<-start.index[-length(start.index)]
    }
    remove.index<-NULL
    for (i in 1:(n-1)){
      next.word<-start.index+i
      new.matches<-which(text.split[next.word]==ngram.split[i+1])
      if (length(new.matches>0)){  
        start.index<-start.index[new.matches]
        remove.index<-c(remove.index, next.word[new.matches])
      } else {
        start.index<-NULL
      }
    }
    to.replace<-paste(ngram.split, collapse='_')
    if (!is.null(start.index)){
      text.split[start.index]<-to.replace
      text.split<-text.split[-remove.index]
    }
  }
  return.text<-paste(text.split, collapse=' ')
  return(return.text)
}

hardClean<-function(corpus){
  corp.names<-names(corpus)
  text.vector<-NULL
  library(tm)
  corpus<-tm_map(corpus, tolower)
  for (i in 1:length(corpus)){
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
    text.vector<-c(text.vector, curr.text.clean)
  }
  names(text.vector)<-corp.names
  text.vector.source<-VectorSource(text.vector)
  new.corpus<-SimpleCorpus(text.vector.source)
  return (new.corpus)
}

makeYears.h<-function(date.df){
  year.m<-NULL
  for (i in 1900:2013){
    if (i %in% date.df$Date){
      num.years.AJPA<-which(date.df$Date == i & date.df$Journal == " American Journal of Physical Anthropology")
      num.years.CA<-which(date.df$Date == i & date.df$Journal == " Current Anthropology,")
      num.years.JFS<-which(date.df$Date == i & date.df$Journal == " Journal of Forensic Sciences")
      to.add<-c(i,length(num.years.AJPA), length(num.years.CA), length(num.years.JFS))
      year.m<-rbind(year.m, to.add)
    }
  }
  colnames(year.m)<-c("Year", "AJPA", "CA", "JFS")
  return (year.m)
}
makeYears.v<-function(date.df){
  year.m<-NULL
  for (i in 1900:2013){
    if (i %in% date.df$Date){
      num.years.AJPA<-which(date.df$Date == i & date.df$Journal == " American Journal of Physical Anthropology")
      to.add.AJPA<-c(i, length(num.years.AJPA), "AJPA")
      num.years.CA<-which(date.df$Date == i & date.df$Journal == " Current Anthropology,")
      to.add.CA<-c(i, length(num.years.CA), "CA")
      num.years.JFS<-which(date.df$Date == i & date.df$Journal == " Journal of Forensic Sciences")
      to.add.JFS<-c(i, length(num.years.JFS), "JFS")
      year.m<-rbind(year.m, to.add.AJPA, to.add.CA, to.add.JFS)
    }
  }
  colnames(year.m)<-c("Year", "Number of Articles", "Journal")
  return (year.m)
}

removeZero<-function(date.df){
  plot.df<-NULL
  for (i in 1:nrow(date.df)){
    if (date.df[i,2] != 0){
      to.add<-date.df[i,]
      plot.df<-rbind(plot.df, to.add)
    }
  }
  return(plot.df)
}

makeCollocateList<-function(ed.corpus, target, horizon){
  collocate.list<-NULL
  #corpus.vector<-unlist(sapply(ed.corpus, '[', "content"))
  for (i in 1:length(ed.corpus)){
    current.text<-ed.corpus[[i]][[1]]
    current.text.list<-strsplit(current.text, " ")
    current.text.list<-unlist(current.text.list)
    if (target %in% current.text.list){
      target.list<-which(current.text.list == target)
      #target.list<-grep(target, current.text.list)
      for (j in 1:length(target.list)){
        start.col<-target.list[j]-horizon
        if (start.col<1){
          start.col<-1
        }
        end.col<-target.list[j]+horizon
        if (end.col>length(current.text.list)){
          end.col<-length(current.text.list)
        }
        collocate.list<-c(collocate.list, current.text.list[c(start.col:(target.list[j]-1),(target.list[j]+1):end.col)])
        #collocate.list<-c(collocate.list, current.text.list[c(start.col:end.col)])
      }
    }
  }
  return(collocate.list)
}
    
makeCollocateTable<-function(collocate.list, min){
  collocate.table<-NULL
  unique.words<-levels(factor(collocate.list))
  for (i in 1:length(unique.words)){
    if (unique.words[i] != ""){
      num.col<-length(which(collocate.list == unique.words[i]))
      collocate.table<-c(collocate.table, num.col)
    }
  }
  bad.words<-which(unique.words=="")
  if(length(bad.words>0)){
    unique.words<-unique.words[-bad.words]
  }
  names(collocate.table)<-unique.words
  below.threshold<-which(collocate.table < min)
  collocate.table<-collocate.table[-below.threshold]
  return(collocate.table)
}


scaleTable<-function(collocate.table){
  numeric.table<-as.numeric(collocate.table)
  scaling<-sum(numeric.table)
  numeric.table<-numeric.table/scaling
  return.table<-rbind(collocate.table, numeric.table)
  rownames(return.table)<-NULL
  return(return.table)
}

getExpected<-function(full.corpus){
  full.dtm<-DocumentTermMatrix(full.corpus, control=list(wordLengths=c(1,Inf)))
  full.matrix<-as.matrix(full.dtm)
  expectedProbs<-colSums(full.matrix)
  scaling<-sum(expectedProbs)
  expectedProbs<-expectedProbs/scaling
  names(expectedProbs)<-colnames(full.matrix)
  return(expectedProbs)
}

getMDWs<-function(obs, freq, total.group.size, alpha){
  MDW.table<-NULL
  words<-names(obs)
  for (j in 1:length(obs)){
    present<-c(obs[j], as.integer(freq[j]*total.group.size))
    absent<-total.group.size-present
    contingency.table<-rbind(present,absent)
    #error control
    word.p<-tryCatch(fisher.test(contingency.table, alternative="g")$p.value, error=function(cond){return(1)})
    if (word.p<alpha){
      obs.exp<-present[1]/(freq[j]*total.group.size)
      curr.MDW<-c(words[j], obs[j], obs.exp, word.p)
      MDW.table<-cbind(MDW.table, curr.MDW)
    } 
  }
  MDW.table.sort<-MDW.table[,order(-as.numeric(MDW.table[3,]))]
  return(MDW.table.sort)
}

#pass a raw Corpus to function, with target list as a vector of words and horizon to get 
#most distinctive words within the horizon of target
#function returns a list containing $by.group (an MDW table containing duplicate values across
#targets) and $across.groups (an MDW table containing MDWs calculated as the aggreate of all collocates of all
#words in the target list)
#new functionality provides for target ngrams of n length which are collapsed in the corpus during analysis

getCollocateMDW<-function(full.corpus, target.list, horizon.value, min=2, alpha=0.05, stopwords=T, collapse.collocates=T, output.filename="TargetCols"){
  library(tm)
  print("Cleaning Corpus")
  full.corpus<-hardClean(full.corpus)
  full.corpus<-tm_map(full.corpus, content_transformer(stripWhitespace))
  all.collocates<-NULL
  individual.collocates<-NULL
  excluded.targets<-NULL
  
  print("Joining Ngrams")
  
  corpus.vector<-lapply(seq(1, length(full.corpus), by=1), function(x) full.corpus[[x]]$content)
  for (i in 1:length(target.list)){
    target.split<-unlist(strsplit(target.list[i], " "))
    if (length(target.split>1)){
      corpus.vector<-lapply(corpus.vector, collapseNgram, ngram=target.list[i])
      target.list[i]<-paste(target.split, collapse='_')
    }
  }
  full.corpus<-SimpleCorpus(VectorSource(corpus.vector))
  remove(corpus.vector)
  print("Calculating Expected Table")
  expected.table<-getExpected(full.corpus)
  
  #calculate MDWs for each target in list
  for (i in 1:length(target.list)){
    print(paste("Now searching collocates of: ", target.list[i], sep=''))
    print("Calculating Collocates")
    collocate.list<-makeCollocateList(full.corpus, target.list[i], horizon.value)
    if (length(collocate.list>0)){
      all.collocates<-c(all.collocates, collocate.list)
      collocate.table<-makeCollocateTable(collocate.list, min)
      total.cols<-sum(collocate.table)
      col.exp.index<-which(names(expected.table) %in% names(collocate.table))
      col.expected<-expected.table[col.exp.index]
      collocate.table.sort<-collocate.table[order(names(collocate.table))]
      expected.table.sort<-col.expected[order(names(col.expected))]
      collocate.MDW<-getMDWs(collocate.table.sort, expected.table.sort, total.cols, alpha)
      #print(collocate.MDW)
      if(!is.null(collocate.MDW)){
        collocate.MDW<-as.matrix(collocate.MDW)
        print(paste("Found", as.numeric(ncol(collocate.MDW)), "distinctive collocates"), sep=' ')
        target.row<-rep(target.list[i], ncol(collocate.MDW))
        collocate.MDW<-rbind(collocate.MDW, target.row)
        obs.values<-as.integer(collocate.MDW[2,])
        obs.exp<-as.integer(collocate.MDW[3,])
        adj.obs<-obs.values/total.cols
        collocate.MDW<-rbind(collocate.MDW, adj.obs)
        #following code creates obs/exp per collocate in group
        #exp.MDW.index<-which(names(col.expected) %in% collocate.MDW[1,])
        #exp.MDW<-col.expected[exp.MDW.index]
        #exp.MDW<-exp.MDW[order(names(exp.MDW))]
        #collocate.MDW.sort<-collocate.MDW[,order(collocate.MDW[1,])]
        #obs.values<-as.integer(collocate.MDW.sort[2,])
        #obs.exp<-obs.values/exp.MDW
      
        individual.collocates<-cbind(individual.collocates, collocate.MDW)
      } else {
        print("Found no distinctive collocates")
        excluded.targets<-c(excluded.targets, target.list[i])
      }
    } else { 
      print("Found no distinctive collocates")
      excluded.targets<-c(excluded.targets, target.list[i])
    }
    write.csv(individual.collocates, file="CollocateTemp.csv")
  }
  words<-individual.collocates[1,]
  nObs<-as.numeric(individual.collocates[2,])
  ObsAdj<-as.numeric(individual.collocates[6,])
  ObsExp<-as.numeric(individual.collocates[3,])
  pValue<-as.numeric(individual.collocates[4,])
  Group<-individual.collocates[5,]
  #build table and output here
  
  if(collapse.collocates){
    #assign each collocate a group based on the Adjusted Obs value, record all groups each collocate belongs to
    obs.raw<-NULL
    obs.adj<-NULL
    p.value<-NULL
    obs.exp<-NULL
    assigned.group<-NULL
    all.groups<-NULL
    unique.collocates<-levels(factor(words))
    for (i in 1:length(unique.collocates)){
      word<-unique.collocates[i]
      word.index<-which(words==word)
      potential.raw.obs<-nObs[word.index]
      potential.obs.exp<-ObsExp[word.index]
      potential.obs<-ObsAdj[word.index]
      potential.p<-pValue[word.index]
      potential.groups<-Group[word.index]
      if (length(word.index) == 1){
        obs.adj<-c(obs.adj, potential.obs)
        obs.raw<-c(obs.raw, potential.raw.obs)
        obs.exp<-c(obs.exp, potential.obs.exp)
        p.value<-c(p.value, potential.p)
        assigned.group<-c(assigned.group, potential.groups)
        all.groups<-c(all.groups, potential.groups)
      } else {
        assign.index<-which(potential.obs == max(potential.obs))
        assign.index<-assign.index[1]
        obs.adj<-c(obs.adj, potential.obs[assign.index])
        obs.exp<-c(obs.exp, potential.obs.exp[assign.index])
        obs.raw<-c(obs.raw, sum(potential.raw.obs))
        p.value<-c(p.value, potential.p[assign.index])
        assigned.group<-c(assigned.group, potential.groups[assign.index])
        all.groups<-c(all.groups, paste(potential.groups, collapse='|'))
      }
    }
    words<-unique.collocates
    nObs<-obs.raw
    obsExp<-obs.exp
    nObsAdj<-obs.adj
    AssignedGroup<-assigned.group
    PotentialGroups<-all.groups
    
    individual.collocates<-data.frame(words, nObs, nObsAdj, obsExp, AssignedGroup, PotentialGroups) 
  } else {
    nObsAdj<-ObsAdj
    individual.collocates<-data.frame(words, nObs, nObsAdj, ObsExp, pValue, Group)
  }
  by.group<-individual.collocates[order(individual.collocates$nObsAdj, decreasing=T),]
  if (stopwords){
    remove.words<-stopwords("English")
    remove.index<-which(by.group$words %in% remove.words)
    by.group<-by.group[-remove.index,]
  }
  
  
  #calculate MDWs for all collocates of all targets
  #depreciated due to utmost confusion
#   collocate.table<-makeCollocateTable(all.collocates, min)
#   total.cols<-sum(collocate.table)
#   col.exp.index<-which(names(expected.table) %in% names(collocate.table))
#   col.expected<-expected.table[col.exp.index]
#   collocate.table.sort<-collocate.table[order(names(collocate.table))]
#   expected.table.sort<-col.expected[order(names(col.expected))]
#   collocate.MDW<-getMDWs(collocate.table.sort, expected.table.sort, total.cols, alpha)
#   words<-collocate.MDW[1,]
#   nObs<-as.integer(collocate.MDW[2,])
#   pValue<-as.numeric(collocate.MDW[3,])
#   unique.collocates<-data.frame(words, nObs, pValue)
#   unique<-unique.collocates[order(unique.collocates$nObs, decreasing=T),]
#   group<-NULL
#   for (i in 1:nrow(unique)){
#     group.index<-which(as.character(by.group[,1]) == as.character(unique[i,1]))
#     if (length(group.index)==0){
#       assign.group<-NA
#     } else {
#       assign.group<-as.character(by.group[group.index,4])
#     }
#     group<-c(group, assign.group)
#   }
#   unique$AssignedGroup<-group
#   across.groups<-unique
#   
#   if (stopwords){
#     remove.words<-stopwords("English")
#     remove.index<-which(across.groups$words %in% remove.words)
#     across.groups<-across.groups[-remove.index,]
#   }
  
  #combine both into a list
  if (is.null(excluded.targets)){
    excluded.targets<-"NULL"
  }
  collocate.combined<-list(by.group, excluded.targets)
  names(collocate.combined)<-c("by.group", "excluded.targets")
  #see above for depreciation of "across.groups"
  #collocate.combined<-list(by.group, across.groups, excluded.targets)
  #names(collocate.combined)<-c("by.group", "across.groups", "excluded.targets")
  write.csv(by.group, file=paste(output.filename, "FinalCollocateTable.csv", sep="_"))
  #return(collocate.combined)
  return(by.group)
} 


stemCollocateList<-function(collocate.table, stemList, groupCol){
  match.list<-NULL
  for (i in 1:length(stemList)){
    curr.target<-unlist(strsplit(stemList[i], ''))
    match.list<-c(match.list, paste(curr.target[1:2], collapse=''))
  }
  stemmedGroups<-NULL
  for (i in 1:nrow(collocate.table)){
    to.match<-unlist(strsplit(as.character(collocate.table[i,groupCol]),''))
    match.index<-which(match.list==paste(to.match[1:2], collapse=''))
    if (length(match.index)!=0){
      stemmedGroups<-c(stemmedGroups, stemList[match.index])
    } else {
      stemmedGroups<-c(stemmedGroups, NA)
    }
  }
  collocate.table$AssignedGroup<-stemmedGroups
  return(collocate.table)
}

collocateNetworks<-function(collocate.table, filename){
  Source<-collocate.table$Group
  Target<-collocate.table$words
  Weight<-collocate.table$ObsExp
  nObsAdj<-collocate.table$nObsAdj
  Type<-rep("Directed", length(Source))
  edge.table<-data.frame(Source, Target, Weight, nObsAdj, Type)
  write.csv(edge.table, filename)
}

#addition for ResourceLit project - create dispersion plots of collocates within texts of a corpus
makeDispersionPlot<-function(target.text, targets, line.color){
  target.text<-unlist(strsplit(target.text, ' '))
  target.table<-data.frame(which(target.text %in% targets), rep(1, length(which(target.text %in% targets))))
  colnames(target.table)<-c("Position", "Height")
  target.plot<-ggplot(target.table, aes(x=Position, y=Height))+geom_vline(color=line.color, aes(xintercept=Position))+ylim(0,1)+xlim(0,length(target.text))+theme_bw()+theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background=element_blank())
  return(target.plot)
}


dispersionPlots<-function(target.corpus, targets, plot.color, corpus.label){
  corpus.clean<-hardClean(target.corpus)
  corpus.texts<-lapply(corpus.clean, function(x) paste(x, collapse=" "))
  corpus.texts<-lapply(corpus.texts, function(x) unlist(strsplit(x, " ")))
  blank.index<-which(corpus.texts=="")
  if(length(blank.index>0)){
    corpus.texts<-corpus.texts[-blank.index]
  }
  target.plots<-lapply(corpus.texts, function(x) makeDispersionPlot(x, targets, plot.color))
  library(gridExtra)
  filename<-paste(corpus.label, "DispersionPlot_Top_50_Perc_TopObsExpTerms.pdf", sep="_")
  all.plots<-marrangeGrob(target.plots, ncol=1, nrow=length(target.plots), top=corpus.label)
  pdf(filename, height=.5*length(target.plots), width=12)
  print(all.plots)
  dev.off()
  return(target.plots)
}

multiDispersionPlots<-function(corpus.list, target.list, color.list, label.list, all.plot.filename){
  all.tables<-mapply(function(x,y,z,a) dispersionPlots(x,y,z,a), corpus.list, target.list, color.list, label.list)
  library(gridExtra)
  plot.grid<-marrangeGrob(all.tables, ncol=2, nrow=(length(all.tables)/2), as.table=FALSE, top="Dispersion Plots")
  pdf(all.plot.filename, width=(15*2), height=(.5*(length(all.tables)/2)))
  print(plot.grid)
  dev.off()
}
