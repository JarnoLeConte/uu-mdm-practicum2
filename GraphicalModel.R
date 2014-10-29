# igraph library
library(igraph)

# model
data <- c(3,17,4,2,176,197,293,23)
margin <- c(2,2,2)
arr <- array(data, margin)
table <- as.table(arr)

# graph
size <- length(dim(table))
adjm <- matrix(0, size, size)

# test search
test <- function() {  
  gm.search(table, adjm, forward=T, backward=T, score="bic") 
} 


gm.search = function(table, adjm, forward, backward, score) {
  
  # XXX hack to explicit cast adjm to be an adjeceny matrix
  # which will not be flatten in the apply function
  adjm <- get.adjacency(graph.adjacency(adjm, mode="undirected"))
  
  # calculate AIC or BIC score of the current model
  current.score <- gm.score(table, adjm, score)
    
  # generate neighbor models. 
  neighbors <- gm.neighbors(adjm, forward, backward)
  
  # calculate AIC or BIC scores of the neighbors
  neighbor.scores <- sapply(neighbors, function(adjm) { 
    list(adjm = adjm, score = gm.score(table, adjm, score))
  });

  # bets neighbor score
  best.neighbor <- neighbor.scores[,which.min(neighbor.scores[2,])]
  
  # if a neighbor has a better score, recurse on it!
  if (best.neighbor$score < current.score) {
    gm.search(table, best.neighbor$adjm, forward, backward, score)
  } else {
    adjm  # DONE! return current model
  }  
}

# generate neighbor models by adding or removing edges to the given graph
gm.neighbors <- function(adjm, forward, backward) {
  
  # determine which edges need to be considered in order to create new neighbor models. 
  # Either add or remove an edge from the graph.
  # i.e. edges are determined by all combinations of 2 vertices 
  considerEdges <- combn(1:nrow(adjm), 2)
  
  # generate neighbor models by (possibly) adding or removing a single edge.
  neighbors <- apply(considerEdges, 2, function(e) { graph.edge.flip(adjm, e, forward, backward) })
  
  # filter the neighbor models which are unmodified, 
  # i.e. models which are the same a the current model 
  neighbors <- Filter(function(m) { !identical(adjm, m) }, neighbors)
  
  # only allow chordal graphs, so remove the models which don't satisfy this constraint
  # i.e. graphs not containing chordless cycles of length > 3
  neighbors <- Filter(graph.is.chordal, neighbors)
  
  neighbors
}

# calculate log linear model and score by AIC or BIC
gm.score <- function(table, adjm, score) { 
  
  # find cliques
  cliques <- graph.find.cliques(c(), 1:nrow(adjm), c(), adjm, list())
  
  # log linear
  model <- loglin(table, cliques, print=T, iter=100)
  
  # score AIC or BIC
  if (score == "aic") gm.score.aic(table, model) else gm.score.bic(table, model)
}

# calculate AIC score of a model
gm.score.aic <- function(table, model) {
  deviance <- model$lrt          # deviance

  maxDF    <- prod(dim(table))
  uterms   <- maxDF - model$df   # number of parameters (u-terms)

  deviance + 2 * uterms          # AIC
}

# calculate BIC score of a model
gm.score.bic <- function(table, model) {
  deviance <- model$lrt          # deviance
  
  maxDF    <- prod(dim(table))
  uterms   <- maxDF - model$df   # number of parameters (u-terms)
  N        <- sum(table)         # number of observations
  
  deviance + log(N) * uterms   # AIC
}



# modify the given edge by adding or removing it from the adjeceny matrix.
# adding the edge is only be done when the forward flag is True
# removing the edge is only be done when the backward flag is True
# if the flag don't matches, the original adjeceny matrix is returned
graph.edge.flip <- function(adjm, edge, forward, backward) {
  
  # check if edge already exists in the graph
  exists <- adjm[edge[[1]],edge[[2]]] == 1
  
  # remove or add the current edge from the adjeceny matrix by flipping the bit
  if ((forward && !exists) || (backward && exists)) {
    adjm[edge[[1]],edge[[2]]] = 1 - exists
  } 
  
  # return a (possibly modified) adjecency matrix
  adjm
}

# check if the graph is chordal
# i.e. check absence of chordless cycles with length > 3
graph.is.chordal <- function(adjm) {
  
  # transform adjeceny matrix to real graph
  graph <- graph.adjacency(adjm, mode="undirected")
  
  # check absence of chordless cycles with length > 3
  is.chordal(graph)$chordal  
}

# given code to find cliques in the graph
graph.find.cliques <- function (R,P,X,graph,cliques) {
  
  # helper function
  neighbors <- function (graph,node) 
  {
    nnodes <- dim(graph)[2]
    (1:nnodes)[graph[node,]==1]
  }
  
  # helper function
  post.process <- function (cliques) 
  {
    unique(lapply(cliques,sort))
  }
  
  if (length(P)==0 & length(X)==0) {cliques <- list(R)}
  else {
    pivot <- P[sample(length(P),1)]
    for(i in 1:length(P)){
      pivot.nb <- neighbors(graph,pivot)
      if(!is.element(P[i],pivot.nb)){
        P.temp <- setdiff(P,P[i])
        R.new <- union(R,P[i])
        P.new <- intersect(P.temp,neighbors(graph,P[i]))
        X.new <- intersect(X,neighbors(graph,P[i]))
        cliques <- c(cliques,find.cliques(R.new,P.new,X.new,graph,list()))
        X <- union(X,P[i])}
    }}
  cliques
}



