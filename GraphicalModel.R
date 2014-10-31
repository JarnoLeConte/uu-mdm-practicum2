

# call gm.search nstart times with different initial matrices
# this will give reliable results instead of calling gm.search once.
gm.restart = function(nstart, prob, seed, table, forward, backward, score) {
  
  # col/row size
  size <- length(dim(table))
  
  # generate symmetric matrix with initial edges given by the probability 'prob'
  matrix.rand <- function(i) {
    set.seed(seed + i)
    elms <- sample(c(0,1), size=size*size, replace=T, prob=c(1-prob, prob))
    m <- matrix(elms, size, size)
    ind <- lower.tri(m) 
    m[ind] <- t(m)[ind]
    m
  }
  
  # single call to gm.search by passing a new generated graph, which is represented 
  # by an adjeceny matrix representing the initial edges
  search <- function(i) {
    adjm <- matrix.rand(i)
    gm.search(table, adjm, forward, backward, score)
  }
  
  # multiple calls to gm.search, each time with a new generated matrix
  searches <- sapply(1:nstart, search)
  
  # find best search
  best.search <- searches[,which.min(searches[2,])]
  
  # return best search
  list(
    search = best.search, 
    call = match.call()
  )
}

# perform a single search given a graph which is represented by an adjecency matrix
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
  
  # best neighbor score
  if (length(neighbors) > 0)
    best.neighbor <- neighbor.scores[,which.min(neighbor.scores[2,])]
  
  # recurse if a neighbor has a better score or stop otherwise
  if (length(neighbors) > 0 && best.neighbor$score < current.score) {
    
    # some neighbor has a better score, recurse on it!
    res <- gm.search(table, best.neighbor$adjm, forward, backward, score)
    
    # append step description to trace.
    diff <- best.neighbor$adjm - adjm;
    diff.indices <- which(as.matrix(diff)==1, arr.ind = T)
    message <- character();
    if (sum(diff) > 0) {
      message <- "Added: "    
    } else {
      message <- "Removed: "
    }    
    
    message <- paste(message, diff.indices[1], " - ", diff.indices[2], " (score= ", best.neighbor$score, ")")
    res$trace <- rbind(message, res$trace)
    
    res
    
  } else {
      
    # DONE! return current model
    list(
      model = graph.find.cliques(adjm),
      score = current.score,
      trace = data.frame(test=character(), stringsAsFactors=FALSE),
      call = match.call()
    )
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
  cliques <- graph.find.cliques(adjm)
  
  # log linear
  model <- loglin(table, cliques, print=F, iter=100)
  
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


#########
# GRAPH #
#########

# igraph library
library(igraph)

# modify the given edge by adding or removing it from the adjeceny matrix.
# adding the edge is only be done when the forward flag is True
# removing the edge is only be done when the backward flag is True
# if the flag don't matches, the original adjeceny matrix is returned
graph.edge.flip <- function(adjm, edge, forward, backward) {
  
  # check if edge already exists in the graph
  exists <- adjm[edge[[1]],edge[[2]]] == 1
  
  # remove or add the current edge from the adjeceny matrix by flipping the bit
  # because the matrix is symmetric, we must do this on 2 positions.
  if ((forward && !exists) || (backward && exists)) {
    adjm[edge[[1]],edge[[2]]] = 1 - exists
    adjm[edge[[2]],edge[[1]]] = 1 - exists
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

graph.find.cliques = function(adjm) {
  
  # transform adjeceny matrix to real graph
  graph <- graph.adjacency(adjm, mode="undirected")
  
  # cliques in the graph
  maximal.cliques(graph)
}


########
# TEST #
########

testCare <- function() {
  
  # observed data (care, survival, clinic)
  data <- c(3,17,4,2,176,197,293,23)
  a <- array(data, c(2,2,2))
  
  # start search
  gm.restart(nstart=5, prob=0.5, seed=2, as.table(a), forward=T, backward=T, score="aic") 
}

testCoronary <- function() {
  
  #src <- "/edu/mdm/mdm2/uu-mdm-practicum2/coronary.txt" # Mathijs
  src <- "/Users/Jarno/dev/school/mdm/uu-mdm-practicum2/coronary.txt" # Jarno
  
  # observed data
  coronary.dat <- read.table(src, header=T)
  
  # start search
  gm.restart(nstart=5, prob=0.5, seed=2, table(coronary.dat), forward=T, backward=T, score="bic") 
}

testRHC <- function() {
  
  #src <- "/edu/mdm/mdm2/uu-mdm-practicum2/rhc-small.txt" # Mathijs
  src <- "/Users/Jarno/dev/school/mdm/uu-mdm-practicum2/rhc-small.txt" # Jarno
  
  # observed data
  rhc.dat <- read.csv(src, header=T)
  
  # start search
  gm.restart(nstart=1, prob=0, seed=2, table(rhc.dat), forward=T, backward=T, score="bic") 
}

