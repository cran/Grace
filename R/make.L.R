# This function makes graph Laplacian matrices from adjacency matrices
make.L <- function(adj, normalize.Laplacian = FALSE){
  if(isSymmetric(adj) == FALSE){
    stop("Error: The adjacency matrix needs to be symmetric.")
  }
  L <- -adj
  diag(L) <- 0
  diag(L) <- -rowSums(L)
  if(normalize.Laplacian){
    diag(L)[diag(L) == 0] <- 1
    L <- diag(1 / sqrt(diag(L))) %*% L %*% diag(1 / sqrt(diag(L)))  # Normalize L
  }
  return(L)
}
