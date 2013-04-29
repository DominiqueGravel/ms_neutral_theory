########################################################################
# Functions to generate networks ("graphs") for spatial ecology.
#
# Stores the graph in a adjacency matrix
#
# Mains functions included:
# - geograph: returns a random geometric graph
# - geotree: returns a random geometric tree
# - lattice: returns a lattice type of spatial graph
# - connected: returns a fully connected graph
# - plot_spatial: custom function to plot the graphs
#
# Associated functions:
# - isConnected: tests if the graph is fully connected
# - testCon: function used in isConnected
# - spt: generates a graph of shortest path given an input graph 
#
# By: Dominique Gravel (dominique_gravel@uqar.ca) & Philippe Desjardins-Proulx (philippe.d.proulx@gmail.com)
# April 2013
########################################################################

geograph_fn = function(n = 64, r = 0.32) {
  # Generates a random geometric graph.
  #
  # Args:
  #   n: Number of vertices.
  #   r: Threshold distance to connect vertices.
  #
  # Returns:
  #   A list with the xy coordinates and the adjacency matrix.
  # Note: the algorithm tests if the graph is connected
  attempts = 0
  repeat {
    attempts <- attempts + 1
	xy = cbind(runif(n), runif(n))
    distMat = as.matrix(dist(xy, method = 'euclidean', upper = T, diag = T))
    adjMat = matrix(0, nr = n, nc = n)
    adjMat[distMat < r] = 1
    diag(adjMat) = 0
    if (isConnected(adjMat)) {
      return(list(xy, adjMat))
    }
  }
}

##############################
geotree_fn = function(n = 64, r = 0.32) {
  # Generates a random geometric tree.
  #
  # Args:
  #   n: Number of vertices.
  #   r: Threshold distance to connect vertices.
  #
  # Returns:
  #   A list with the xy coordinatesa nd the adjacency matrix
  source <- sample(1:n, 1)
  g <- geograph_fn(n, r)
  return(list(g[[1]], spt(g[[2]], source)))
}

##############################
lattice_fn = function(n) {
	# Generates a lattice type of graph
 	#
 	# Args:
 	#   n: number of cells in the lattice
 	#
 	# Returns:
 	#   A list with the xy coordinates and the adjacency matrix  
	X = seq(0,1,by = 1/(n^0.5-1))
	Y = seq(0,1,by = 1/(n^0.5-1))
	XY = expand.grid(X,Y)
	distMat = as.matrix(dist(XY,method = "euclidean", upper = T, diag = T))
	adjMat = matrix(0, nr=n, nc=n)
	adjMat[distMat <= 1/(n^0.5-1)*(1+1e-10)] = 1
	diag(adjMat) = 0
	return(list(XY,adjMat)) 	
}

##############################
connected_fn = function(n) {
	# Generates a connected graph
 	#
 	# Args:
 	#   n: number of cells in the graph
 	#
 	# Returns:
 	#   A list with the xy coordinates and the adjacency matrix  
	X = sin(c(0:(n-1))*2*pi/(n-1))
	Y = cos(c(0:(n-1))*2*pi/(n-1))
	adjMat = matrix(1,nr = n, nc = n)
	diag(adjMat) = 0
	XY = cbind(X,Y)
	return(list(XY,adjMat)) 	
}

##############################
isConnected = function(adjMat) {
  # Tests if the graph is fully connected (i.e.: there is a path between all nodes).
  #
  # Args:
  #   adjMat: An adjacency matrix (made of boolean values).
  #
  # Returns:
  #   TRUE if the graph is connected, FALSE otherwise.
  diag(adjMat) = 0
  n = nrow(adjMat)
  for (v in 1:n) {
    inPath = vector('logical', n)
    inPath[v] = T
    inPath = testCon(adjMat, inPath, v)
    if (sum(inPath == T) < n) {
      return(F)
    }
  }
  return(T)
}

##############################
testCon = function(adjMat, inPath, v) {
  # Helper recursive function for 'isConnected'.
  for (i in 1:nrow(adjMat)) {
    if (adjMat[v, i] && inPath[i] == F) {
      inPath[i] = T
      inPath = testCon(adjMat, inPath, i)
    }
  }
  return(inPath)
}

##############################
spt = function(adjMat0, source) {
  # Generates the graph of the shortest path tree using the Dijkstra
  # algorithm for a given graph and source vertex.
  #
  # Args:
  #   adjMat0: The adjacency matrix for the graph.
  #   source: The starting vertex for the algorithm.
  #
  # Returns:
  #   The shortest path tree as an adjacency matrix (a matrix of bool).
  n = nrow(adjMat0)
  distance = rep(9999, n)
  distance[source] = 0
  previous = rep(0, n)
  q = 1:n

  while (length(q) != 0) {
    u = q[1]
    for (i in q) {
      if (distance[i] < distance[u]) {
        u = i
      }
    }
    q = setdiff(q, u)
    for (i in q) {
      x = distance[u]
      if (adjMat0[u, i]) {
        x = x + 1
      } else {
        x = x + 9999
      }
      if (x < distance[i]) {
        distance[i] = x
        previous[i] = u
      }
    }
  }
  
  # Build the adjacency matrix:
  adjMat = matrix(0, nr = n, nc = n)
  for (i in 1:n) {
    adjMat[i, previous[i]] <- 1
    adjMat[previous[i], i] <- 1
  }
  return(adjMat)
}

##############################
spm = function(adjMat) {
  # Generates the shortest path distance matrix for unweighted graphs using
  # the Dijkstra algorithm.
  #
  # Args:
  #   adjMat: The adjacency matrix for the unweighted graph (should be
  #           filled with 0s and 1s).
  #
  # Returns:
  #   A matrix with the distance between all pair of vertices.
  m = as.matrix(spd(adjMat, 1))
  for (i in 2:nrow(adjMat)) m = cbind(m, spd(adjMat, i))
  return(m)
}

##############################
spd = function(adjMat, source) {
  # Generates the shortest path distances starting with a given source.
  #
  # Args:
  #   adjMat: The adjacency matrix for the graph.
  #   source: The starting vertex for the algorithm.
  #
  # Returns:
  #   A vector with the shortest path distance.
  n = nrow(adjMat)
  distance = rep(Inf, n)
  distance[source] = 0
  previous = rep(0, n)
  q = 1:n

  while (length(q) != 0) {
    u = q[1]
    for (i in q) {
      if (distance[i] < distance[u]) {
        u = i
      }
    }
    q = setdiff(q, u)
    for (i in q) {
      x = distance[u]
      if (adjMat[u, i] == 1) {
        x = x + 1
      } else {
        x = Inf
      }
      if (x < distance[i]) {
        distance[i] = x
        previous[i] = u
      }
    }
  }
  return(distance)
}

##############################

plot_spatial = function(spatial_graph, vec.col) {
  	# Plots a spatial graph 
  	#
  	# Args:
  	#   spatia_graph: the output of one of the spatial graphs
  	#
	x11(height = 5.5, width = 6)
	par(mar=c(5,6,2,1))
	XY = spatial_graph[[1]]
	adjMat = spatial_graph[[2]]

	plot(XY[,1],XY[,2],xlab = "X", ylab = "Y",cex.lab = 1.5, cex.axis = 1.25)
	adjVec = stack(as.data.frame(adjMat))[,1]
	XX = expand.grid(XY[,1],XY[,1])
	YY = expand.grid(XY[,2],XY[,2])
	XX = subset(XX,adjVec==1)
	YY = subset(YY,adjVec==1)
	arrows(x0 = XX[,1],x1=XX[,2],y0 = YY[,1], y1 = YY[,2], length = 0,lwd = 0.1, col = "grey")
	points(XY[,1],XY[,2],pch=21,bg=vec.col)
}


##############################
# Examples

#plot_spatial(connected(n = 25))
#plot_spatial(lattice(n = 25))
#plot_spatial(geograph(n = 25, r = 0.3))
#plot_spatial(geotree(n = 25, r = 0.3))



