
# graph represented by an adjency matrix
adjmatrix <- matrix(0,6,6)





gm.search = function(table, graph, forward, backward, score) {
  # XXX
  
  # return
#   list(
#     $model = list(c(1,2,3), c(1,2,3)), 
#     $score = "bic", 
#     $trace = data.frame(), 
#     $call = match.call()
#   )
}

iteration <- function() {
  
  
}

model.aic <- function(model) {
  dev <- model$lrt              # deviance
  dim <- length(dim(model$fit)) # number of parameters (u-terms)
  dev + 2*dim                   # AIC
}




# igraph library

library(igraph)

# create initial graph from adjecency matrix
graph <- graph.adjacency(adjmatrix, mode="undirected")

# add edges (cycle + diagonal)
graph <- graph + edge(1,2)
graph <- graph + edge(2,3)
graph <- graph + edge(3,4)
graph <- graph + edge(4,1)
#graph <- graph + edge(2,4)

# check absence of chordless cycles with length > 3
is.chordal(graph)$chordal  

adjm <- get.adjacency(graph, type="both")

# determine which edges need to be considered in order to create new neighbor models. 
# Either add or remove an edge from the graph.
# i.e. edges are determined by all combinations of 2 vertices 
considerEdges <- combn(1:nrow(adjm), 2)

# generate all possible neighbor models by adding or removing a single edge.
neighbors <- apply(considerEdges, 2, function(e) {
  
  # remove or add the current edge from the adjeceny matrix by flipping the bit.
  adjm[e[[1]],e[[2]]] = 1 - adjm[e[[1]],e[[2]]]
  
  # return a new adjecency matrix containing the modified edge
  adjm
})

# only allow chordal graphs
# i.e. graphs not containing chordless cycles of length > 3
neighbors <- Filter(function(adjm) { 
  graph <- graph.adjacency(adjm, mode="undirected"); 
  is.chordal(graph)$chordal
}, neighbors)

# calculate AIC of neighbors
aics <- sapply(neighbors, function(model) {
  
});
