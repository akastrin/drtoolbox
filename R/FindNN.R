FindNN <-
function(x, k) {
  k <- k + 1
  n <- dim(x)[1]
  D <- as.matrix(dist(x))
  ind <- t(apply(D, 1, order))
  a <- repmat(t(t(1:n)), 1,n-k)
  flat <- a + n * ind[, (k+1):ncol(ind)] - n
  D[as.vector(flat)] <- 0
  D[seq(1, n * n, by = n + 1)] <- 1e-7
  ni <- ind[, 1:k]
  ni <- ni[, 2:ncol(ni)]
  return(list(D, ni))
}
