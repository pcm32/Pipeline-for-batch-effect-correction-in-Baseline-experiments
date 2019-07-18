library(magrittr)
library(stringr)
library(purrr)
library(ggplot2)
library(igraph)
library(gtools)

library(sva)
library(RUVSeq)
library(batchelor)

list_experiments<-function(directory){
  experiments<-list()
  for(file in dir(directory)){
    experiments[[file %>% str_split("-") %>% extract2(1) %>% extract(2:3) %>% paste(collapse="")]] <- get(load(paste0(directory,'/',file)))$rnaseq
  }
  return(experiments)
}

remove_isolated_experiments<-function(experiments,factor){
  batch<-experiments %>% imap(~.y %>% rep(dim(.x)[2])) %>% unlist(use.names=FALSE)
  group<-experiments %>% map(~.[[factor]]) %>% unlist(use.names=FALSE)
  # batch<-NULL; group<-NULL
  # for(i in experiments %>% seq_along){
  #   batch%<>%c(names(experiments)[[i]] %>% rep(dim(experiments[[i]])[2]))
  #   group%<>%c(experiments[[i]][[factor]])
  # }; batch%<>%factor
  groups <- group %>% split(batch)
  
  intersections<-NULL
  for(i in groups){
    for(j in groups){
      intersections%<>%c(length(intersect(i,j)))
    }
  }
  intersections%<>%matrix(length(groups))%<>%set_colnames(names(groups))%<>%set_rownames(names(groups))
  
  intersections %>% graph_from_adjacency_matrix %>% plot
  
  experiments[rownames(intersections)[rowSums(intersections!=0)<=1]]<-NULL
  
  batch<-experiments %>% imap(~.y %>% rep(dim(.x)[2])) %>% unlist(use.names=FALSE)
  group<-experiments %>% map(~.[[factor]]) %>% unlist(use.names=FALSE)
  
  # batch<-NULL; group<-NULL
  # for(i in experiments %>% seq_along){
  #   batch%<>%c(names(experiments)[[i]] %>% rep(dim(experiments[[i]])[2]))
  #   group%<>%c(experiments[[i]][[factor]])
  # }; batch%<>%factor
  
  groups <- group %>% split(batch)
  
  intersections<-NULL
  for(i in groups){
    for(j in groups){
      intersections%<>%c(length(intersect(i,j)))
    }
  }
  intersections%<>%matrix(length(groups))%<>%set_colnames(names(groups))%<>%set_rownames(names(groups))
  
  intersections %>% graph_from_adjacency_matrix %>% plot
  return(experiments)
}

merge_experiments<-function(experiments, log=TRUE, filter=TRUE){
  genes<-experiments %>% map(rownames)
  common.genes<-genes %>% purrr::reduce(intersect)
  filtered.genes<-genes %>% map(setdiff %>% partial(y=common.genes))
  
  data<-experiments %>% map(~.x@assays$data$counts %>% extract(rownames(.x)%in%common.genes,)) %>% purrr::reduce(cbind)
  batch<-experiments %>% imap(~.y %>% rep(ncol(.x))) %>% unlist(use.names=FALSE) %>% factor
  
  # data<-NULL;batch<-NULL
  # for(i in experiments %>% seq_along){
  #   experiments[[i]]->exp
  #   data%<>%cbind(exp %>% assays %>% use_series(counts) %>% extract(rownames(exp)%in%common.genes,))
  #   batch%<>%c(names(experiments)[[i]] %>% rep(dim(exp)[2]))
  # }; batch%<>%factor
  if(filter){
    filter <- data %>% t %>% data.frame %>% split(batch) %>% map(~colSums(.)!=0) %>% purrr::reduce(`&`)
    data%<>%extract(filter,)
  }
  if(log) data%<>%log1p
  # data%<>%correct.batch.effect(batch,method,model,log,model.data=if(is.null(model)) NULL else model.frame(model,vars),k=k)
  return(SummarizedExperiment(
    assays=if(log) list(log_counts=data) else list(counts=data),
    colData=experiments %>% map(colData) %>% smartbind(list=.) %>% set_rownames(experiments %>% map(colnames) %>% unlist) %>% cbind(batch),
    metadata=list(batch=batch)
  ))
}

correct_batch_effect<-function(experiment, model, method=c('ComBat','RUV','MNN'), k){
  log<-experiment@assays$data %>% names %>% switch(log_counts=TRUE, counts=FALSE)
  model.data<-model.frame(model, experiment@colData[all.vars(model)])
  return(SummarizedExperiment(
    assays = list(switch(
      method,
      ComBat = ComBat(experiment@assays$data[[1]], experiment$batch, mod=model.matrix(model, data=model.data)),
      RUV = RUVs(experiment@assays$data[[1]], cIdx=seq_len(nrow(experiment@assays$data[[1]])), k=k, 
                 scIdx=model.data %>% expand.grid %>% apply(1,paste) %>% makeGroups, isLog=log)$normalizedCounts,
      MNN = mnnCorrect(experiment@assays$data[[1]], batch=experiment$batch, k=k)@assays$data$corrected
    )) %>% set_names(if(log) 'corrected_log_counts' else 'corrected_counts'),
    colData = experiment@colData,
    metadata = experiment@metadata
  ))
}