#
# Data Mining 2014
# Assignment 2: Graphical Modeling
#
# Jarno Le Conte (3725154)
# Mathijs Baaijens (3542068)
#

# Learns an undirected graphical model from the input data. Restarts the search multiple
# times with random starting models to prevent getting stuck in local minima.
#
# Args:
#   nstart: Number of times to restart the search.
#   prob: Probability of an edge being present in an initial random model.
#   seed: Random number generator seed.
#   table: Table of observed counts.
#   adjm: Adjacency matrix.
#   forward : Adding edges is allowed.
#   backward : Removing edges is allowed.
#   score : Scoring model.
#
# Returns:
#   Undirected graphical model.
gm.restart = function(nstart, prob, seed, table, forward, backward, score) {
  print("gm.restart");
  
  size <- length(dim(table))
  
  # generate symmetric matrix with initial edges given by the probability 'prob'.
  matrix.rand <- function(i) {
    set.seed(seed + i)
    elms <- sample(c(0,1), size=size*size, replace=T, prob=c(1-prob, prob))
    m <- matrix(elms, size, size)
    ind <- lower.tri(m) 
    m[ind] <- t(m)[ind]  
    m <- m - diag(size);
    m
  }
  
  # single call to gm.search by passing a new generated graph, which is represented 
  # by an adjeceny matrix representing the initial edges.
  search <- function(i) {
    adjm <- matrix.rand(i);
    gm.search(table, adjm, forward, backward, score)
  }
  
  # multiple calls to gm.search, each time with a new generated matrix.
  searches <- sapply(1:nstart, search)
  
  # find best search.
  best.search <- searches[,which.min(searches[2,])]
  
  # return best search.
  list(
    search = best.search, 
    call = match.call()
  )
}

# Learns an undirected graphical model from the input data.
#
# Args:
#   table: Table of observed counts.
#   adjm: Adjacency matrix.
#   forward : Adding edges is allowed.
#   backward : Removing edges is allowed.
#   score : Scoring model.
#
# Returns:
#   Undirected graphical model.
gm.search = function(table, adjm, forward, backward, score) {
  print("gm.search");
  
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
    diff.indices <- which(as.matrix(diff)!=0, arr.ind = T)
    res$trace <- rbind(c(sum(diff) > 0, diff.indices[1], diff.indices[2], best.neighbor$score), res$trace)
    colnames(res$trace) <- c("forward", "start", "end", "score")
    
    res
    
  } else {
      
    # DONE! return current model
    list(
      model = post.process(graph.find.cliques(adjm)),
      score = current.score,
      trace = NULL, 
      call = match.call()
    )
  }  
}

# Generate neighbor models by adding or removing edges to the given graph
#
# Args:
#   adjm: Adjacency matrix.
#   forward : Adding edges is allowed.
#   backward : Removing edges is allowed.
#
# Returns:
#   List of adjeceny matrices
gm.neighbors <- function(adjm, forward, backward) {
  print("gm.neighbors");
  
  # determine which edges need to be considered in order to create new neighbor models. 
  # Either add or remove an edge from the graph.
  # i.e. edges are determined by all combinations of 2 vertices 
  considerEdges <- combn(1:nrow(adjm), 2)
  
  # generate neighbor models by (possibly) adding or removing a single edge.
  neighbors <- apply(considerEdges, 2, function(e) { graph.edge.flip(adjm, e, forward, backward) })
  
  # filter the neighbor models which are unmodified, 
  # i.e. models which are the same a the current model 
  neighbors <- Filter(function(m) { !identical(adjm, m) }, neighbors)
  
  neighbors
}

# Calculate log linear model and score by AIC or BIC
#
# Args:
#   table: Table of observed counts.
#   adjm: Adjacency matrix.
#   score : AIC or BIC
#
# Returns:
#   AIC or BIC score
gm.score <- function(table, adjm, score) { 
  print("gm.score");
  
  # find cliques
  cliques <- graph.find.cliques(adjm)
  
  # log linear
  model <- loglin(table, cliques, print=F, iter=100)
  
  # score AIC or BIC
  if (score == "aic") 
    gm.score.aic(table, model) 
  else 
    gm.score.bic(table, model)
}

# Calculate AIC score of a model
#
# Args:
#   table: Table of observed counts.
#   model: Results from loglin
#
# Returns:
#   AIC score
gm.score.aic <- function(table, model) {
  deviance <- model$lrt          # deviance

  maxDF    <- prod(dim(table))
  uterms   <- maxDF - model$df   # number of parameters (u-terms)

  deviance + 2 * uterms          # AIC
}

# Calculate BIC score of a model
#
# Args:
#   table: Table of observed counts.
#   model: Results from loglin
#
# Returns:
#   BIC score
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

# Modify the given edge by adding or removing it from the adjeceny matrix.
# adding the edge is only be done when the forward flag is True
# removing the edge is only be done when the backward flag is True
# if the flag don't matches, the original adjeceny matrix is returned
#
# Args:
#   adjm: Adjacency matrix.
#   edge: c(i,j) where i and j are column and row indicies
#   forward : Adding edges is allowed.
#   backward : Removing edges is allowed.
#
# Returns:
#   Adjeceny matrix
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


# Find cliques in the graph
#
# Args:
#   adjm: Adjacency matrix.
#
# Returns:
#   List of cliques
graph.find.cliques = function(adjm) {
  size <- dim(adjm)[1]
  post.process(find.cliques(c(), (1:size), c(), adjm, list()))
}

find.cliques <- function (R,P,X,graph,cliques) 
{ 
  neighbors <- function (graph,node) 
  {
    nnodes <- dim(graph)[2]
    (1:nnodes)[graph[node,]==1]
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

post.process <- function (cliques) 
{
  unique(lapply(cliques,sort))
}




########
# TEST #
########

# observed data (care, survival, clinic)
data <- c(3,17,4,2,176,197,293,23)
a <- array(data, c(2,2,2))

#src1 <- "/edu/mdm/mdm2/uu-mdm-practicum2/data/coronary.txt" # Mathijs
src1 <- "/Users/Jarno/dev/school/mdm/uu-mdm-practicum2/data/coronary.txt" # Jarno

# observed data
coronary.dat <- read.table(src1, header=T)

#src2 <- "/edu/mdm/mdm2/uu-mdm-practicum2/data/rhc-small.txt" # Mathijs
src2 <- "/Users/Jarno/dev/school/mdm/uu-mdm-practicum2/data/rhc-small.txt" # Jarno

# observed data
rhc.dat <- read.csv(src2, header=T)


testCare <- function() {
  gm.restart(nstart=5, prob=0.5, seed=2, as.table(a), forward=T, backward=T, score="aic") 
}

testCoronary <- function() {  
  gm.restart(nstart=5, prob=0.5, seed=2, table(coronary.dat), forward=T, backward=T, score="bic") 
}

testRHC <- function() {
  gm.restart(nstart=1, prob=0, seed=2, table(rhc.dat), forward=T, backward=T, score="bic") 
}

testRHC2 <- function() {
  gm.restart(nstart=1, prob=1, seed=2, table(rhc.dat), forward=T, backward=T, score="bic") 
}

testRHC3 <- function() {
  gm.restart(nstart=1, prob=0, seed=2, table(rhc.dat), forward=T, backward=T, score="aic") 
}
testRHC3_B <- function() {
  gm.restart(nstart=1, prob=1, seed=2, table(rhc.dat), forward=T, backward=T, score="aic") 
}

testA1 <- function() {
  gm.restart(nstart=1, prob=0.25, seed=2, table(rhc.dat), forward=T, backward=T, score="aic") 
}
testA2 <- function() {
  gm.restart(nstart=1, prob=0.5, seed=2, table(rhc.dat), forward=T, backward=T, score="aic") 
}
testA3 <- function() {
  gm.restart(nstart=1, prob=0.75, seed=2, table(rhc.dat), forward=T, backward=T, score="aic") 
}

testRestartA1 <- function() {
  gm.restart(nstart=3, prob=0.25, seed=2, table(rhc.dat), forward=T, backward=T, score="aic") 
}
testRestartA2 <- function() {
  gm.restart(nstart=3, prob=0.5, seed=2, table(rhc.dat), forward=T, backward=T, score="aic") 
}
testRestartA3 <- function() {
  gm.restart(nstart=3, prob=0.75, seed=2, table(rhc.dat), forward=T, backward=T, score="aic") 
}

testB1 <- function() {
  gm.restart(nstart=1, prob=0.25, seed=2, table(rhc.dat), forward=T, backward=T, score="bic") 
}
testB2 <- function() {
  gm.restart(nstart=1, prob=0.5, seed=2, table(rhc.dat), forward=T, backward=T, score="bic") 
}
testB3 <- function() {
  gm.restart(nstart=1, prob=0.75, seed=2, table(rhc.dat), forward=T, backward=T, score="bic") 
}

testRestartB1 <- function() {
  gm.restart(nstart=3, prob=0.25, seed=2, table(rhc.dat), forward=T, backward=T, score="bic") 
}
testRestartB2 <- function() {
  gm.restart(nstart=3, prob=0.5, seed=2, table(rhc.dat), forward=T, backward=T, score="bic") 
}
testRestartB3 <- function() {
  gm.restart(nstart=3, prob=0.75, seed=2, table(rhc.dat), forward=T, backward=T, score="bic") 
}

testRestartA1_extra <- function() {
  gm.restart(nstart=9, prob=0.25, seed=2, table(rhc.dat), forward=T, backward=T, score="aic") 
}
testRestartA2_extra <- function() {
  gm.restart(nstart=9, prob=0.50, seed=2, table(rhc.dat), forward=T, backward=T, score="aic") 
}
testRestartA3_extra <- function() {
  gm.restart(nstart=9, prob=0.75, seed=2, table(rhc.dat), forward=T, backward=T, score="aic") 
}

testRestartB1_extra <- function() {
  gm.restart(nstart=9, prob=0.25, seed=2, table(rhc.dat), forward=T, backward=T, score="bic") 
}
testRestartB2_extra <- function() {
  gm.restart(nstart=9, prob=0.5, seed=2, table(rhc.dat), forward=T, backward=T, score="bic") 
}
testRestartB3_extra <- function() {
  gm.restart(nstart=9, prob=0.75, seed=2, table(rhc.dat), forward=T, backward=T, score="bic") 
}