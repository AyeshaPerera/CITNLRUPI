library(statmod)
GL_modified_nodes_weights <- function(lower, upper, 
         number.of.nodes){
  # This module computes the "modified"
  # Gauss Legendre (GL) nodes and weights when the 
  # interval of integration is not necessarily 
  # [-1, 1] but is, instead, [lower, upper].
  # 
  # Inputs
  # lower: lower endpoint of interval of integration
  # upper: upper endpoint of interval of integration
  # number.of.nodes: the number of Gauss Legendre 
  #                  nodes
  # 
  # Packages used: statmod
  # 
  # Output
  # A list with two elements:
  # The first element is GL.modified.nodes
  # which is the vector of "modified" nodes.
  # The second element is GL.modified.weights
  # which is the vector of "modified" weights.
  # 
  # Written by P.Kabaila in Oct 2022
  
  if (upper <= lower){
    stop("upper is less than or equal to lower")
  }
  
  GL.list <- gauss.quad(number.of.nodes, kind = "legendre")
  GL.nodes <- GL.list$nodes
  GL.weights <- GL.list$weights
  
  middle <- (lower + upper) / 2
  half.width <- (upper - lower) / 2
  
  GL.modified.nodes <- middle + half.width * GL.nodes
  GL.modified.weights <- half.width * GL.weights
  
  list(GL.modified.nodes = GL.modified.nodes,
       GL.modified.weights = GL.modified.weights)
  
}