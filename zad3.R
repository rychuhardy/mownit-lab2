library(igraph)
library(cape)
library(compare)

main <-function(filename, s, t, E)
{
  
  M <- read_matrix(filename) # Read matrix from file (Default format is: v1 v2 weight)
  
  G <- graph_from_adjacency_matrix(M, mode="undirected", weighted="weigh") # transform to igraph
  
  G <- set_edge_attr(G, "sem", index=get.edge.ids(G, c(s,t), directed=FALSE), value=E) # set sem - used for drawing
  
  res <- put_sem(G, s, t, E) # Calculates A and B (from equation A*I=B)
  
  A <- matrix(unlist(res[1]), ncol=sqrt(length(unlist(res[1]))), nrow=sqrt(length(unlist(res[1])))) # Extract A matrix
  
  B <- cbind(unlist(res[2])) # Extract B matrix
  
  eg <- solve(A)%*%B # I is A^-1 * B
  
  idx <- 1
  for(i in E(G)) { # Set current to each edge
    G <- set_edge_attr(G, "current", index=i, value=eg[idx])
    idx <- idx+1
  }
  
  plot(G,vertex.shape="circle", layout=layout.circle, edge.label=E(G)$current) # Plot graph with calculated current
  
}

read_matrix <- function(filename) # returns adjacency matrix
{
  table <- read.csv(filename, sep=" ", comment.char = "#", header=FALSE);
  
  max_vertex <- max(table[,1:2])
  
  A <- matrix(c(rep(0,max_vertex*max_vertex)),ncol=max_vertex, nrow=max_vertex) # Create matrix for storing data from file
  
  for(i in 1:nrow(table)) {
    A[table[i,1],table[i,2]] <- table[i,3]
    A[table[i,2],table[i,1]] <- table[i,3]
  }
  
  return(A)
}

remove_repetitions <- function(L) # Helper function since unique() doesn't work with igraphs. Removes repeating cycles from L
{
  for(i in 1:length(L)) {
    j <- i+1
    while(j <= length(L)) {
      if(isTRUE(compare(L[i], L[j], allowAll = TRUE, shorten = FALSE))) {
        L <- L[-j]
        j<- j-1
      }
      j<-j+1
    }
  }
  return (L)
}

find.cycles <- function(graph, k) { # Helper function. Returns all cycles of lenth k (with repetitions)
  ring <- graph.ring(k)
  graph.get.subisomorphisms.vf2(graph, ring)
}


find_cycles <- function(G, i, j) { # Returns list of cycles of length [i,j] in graph G
  
  res <- c()
  for(k in i:j) {
    tmp <- find.cycles(G, k)
    tmp <- remove_repetitions(tmp)
    res <- append(res, tmp)
  }
  return(res)
}

put_sem <-function(G,s,t, E)
{
  cls <- find_cycles(G,3,length(E(G))) # Cycles are used in the KVL
  
  Sys <- matrix(rep(0, length(E(G))*(length(E(G)))), nrow=length(E(G)), ncol=(length(E(G)))) # Stores parameters A (from equation A*I=B)
  cols <- cbind(rep(0, length(E(G)))) # Stores parameters B (from equation A*I=B)
  
  # Values of source
  cols[s] <- E 
  cols[t] <- -E
  
  #Assume directions of current in graph: in edge e(v1,v2) the current goes from min(v1,v2) to max(v1,v2)
  
  #Get equations from KCL
  for (node in V(G)) {
    for(nei in neighbors(G,node)) {
      if(node < nei) {
      Sys[node, get.edge.ids(G, c(node,nei), directed=FALSE)] <- -1 # Node is smaller than nei so the current goes from node to nei
      }
      else {
        Sys[node, get.edge.ids(G, c(node,nei), directed=FALSE)] <- 1 # Opposite of the above
      }
    }
  }
  
  # Remove linearly dependent equations
  Sys <- find_lin(Sys)
  idx <- nrow(Sys)+1
  # The above function removed rows, so restoring the missing rows with zeros
  tmp <- matrix(rep(0, length(E(G))*(length(E(G))-nrow(Sys))), nrow=length(E(G))-nrow(Sys), ncol=(length(E(G))))
  Sys <- rbind(Sys, tmp)
  
  #Get equations from KVL
  for(i in 1:length(cls)) {
    if(idx <= length(E(G))) { # Get only as much equations as needed
    d <- unlist(cls[i]) 
    for(j in 1:length(d)) {
      if(j==length(d)) { # Modulo doesn't work here because indexing starts at 1
        if(d[j] < d[1]) { # The case with edge (last,first) on the cycle list
          Sys[idx,get.edge.ids(G, c(d[j],d[1]), directed=FALSE)] <- get.edge.attribute(G, "weigh", index=get.edge.ids(G, c(d[j],d[1]), directed=FALSE))
        }
        else { 
          Sys[idx,get.edge.ids(G, c(d[j],d[1]), directed=FALSE)] <- get.edge.attribute(G, "weigh", index=get.edge.ids(G, c(d[j],d[1]), directed=FALSE))
        }
      }
      else {
        if(d[j]<d[j+1]) { # the current goes from d[j] to d[j+1]
          Sys[idx,get.edge.ids(G, c(d[j],d[j+1]), directed=FALSE)] <- get.edge.attribute(G, "weigh", index=get.edge.ids(G, c(d[j],d[j+1]), directed=FALSE))
        }
        else {
          Sys[idx,get.edge.ids(G, c(d[j],d[j+1]), directed=FALSE)] <- -get.edge.attribute(G, "weigh", index=get.edge.ids(G, c(d[j], d[j+1]), directed=FALSE))
        }
      }
    }
    }
    idx <- idx+1
  }
  
  return (list(Sys,cols))
}

find_lin<-function(A) # Removes linearly dependent equations (from rows)
  # How this function works is really interesting, unfortunately I haven't worked it out yet :-)
{
  A <- t(A)
  
  q <- qr(A)
  mmat <- A[,q$pivot[seq(q$rank)]]
  
  A<-t(mmat)
  
  return(A)
  
}