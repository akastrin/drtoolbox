\name{IntrinsicDimension}
\alias{IntrinsicDimension}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
IntrinsicDimension(X, method = c("CorrDim", "NearNbDim", "PackingNumbers", "GMST", "EigValue", "MLE"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, method = c("CorrDim", "NearNbDim", "PackingNumbers", 
    "GMST", "EigValue", "MLE")) 
{
    X <- unique(X)
    X <- scale(X)
    method <- match.arg(method)
    if (method == "CorrDim") {
        n <- dim(X)[1]
        D <- FindNN(x = X, k = 5)[[1]]
        r1 <- median(D[D != 0])
        r2 <- max(D[D != 0])
        s1 <- 0
        s2 <- 0
        X <- t(X)
        XX <- apply(X^2, 2, sum)
        onez <- rep(1, n)
        for (i in 1:n) {
            if (i < n) {
                p <- t(X[, i])
                xx <- sum(XX[i])
                xX <- p \%*\% X
                dist <- xx * onez + XX - 2 * xX
                dist <- sqrt(dist[(i + 1):n])
                s1 <- s1 + length(dist[dist < r1])
                s2 <- s2 + length(dist[dist < r2])
            }
        }
        Cr1 <- (2/(n * (n - 1))) * s1
        Cr2 <- (2/(n * (n - 1))) * s2
        no.dims <- (log(Cr2) - log(Cr1))/(log(r2) - log(r1))
    }
    if (method == "NearNbDim") {
        k1 <- 6
        k2 <- 12
        res <- FindNN(X, k2)
        D <- res[[1]]
        ind <- res[[2]]
        Tk <- rep(0, k2 - k1)
        for (k in k1:k2) {
            r <- 1:dim(X)[1]
            c <- ind[, k]
            bla <- cbind(r, c)
            Tk[k - k1 + 1] <- sum(D[bla])
        }
        Tk <- Tk/dim(X)[1]
        no.dims <- (log(Tk[length(Tk)]) - log(Tk[1]))/(log(k2) - 
            log(k1))
    }
    if (method == "PackingNumbers") {
        r <- NULL
        r[1] = 0.1
        r[2] = 0.5
        epsilon <- 0.01
        max.iter <- 20
        done <- 0
        l <- 0
        a <- NULL
        while (done == 0) {
            l <- l + 1
            perm <- sample(dim(X)[1])
            for (k in 1:2) {
                C <- NULL
                for (i in 1:dim(X)[1]) {
                  for (j in 1:length(C)) {
                    if (sqrt(sum((X[perm[i], ] - X[C[j], ])^2)) < 
                      r[k]) {
                      j <- dim(X)[1] + 1
                      break
                    }
                  }
                  if (length(C) == 0 || j < dim(X)[1] + 1) {
                    C <- rbind(C, perm[i])
                  }
                }
                a = c(a, log(length(C)))
            }
            A <- matrix(a, 2)
            no.dims <- -((mean(A[2, ]) - mean(A[1, ]))/(log(r[2]) - 
                log(r[1])))
            if (l > 10) {
                if (1.65 * (sqrt(var(A[1, ])^2 + var(A[2, ])^2)/(sqrt(l) * 
                  log(r[2]) - log(r[1]))) < no.dims * ((1 - epsilon)/2)) {
                  done <- 1
                }
            }
            if (l > max.iter) {
                done <- 1
            }
        }
    }
    if (method == "GMST") {
        gamma <- 1
        M <- 1
        N <- 10
        samp.points <- (dim(X)[1] - 10):(dim(X)[1] - 1)
        k <- 6
        Q <- length(samp.points)
        knnlenavg.vec <- matrix(0, M, Q)
        knnlenstd.vec <- matrix(0, M, Q)
        dvec <- rep(0, 1)
        D <- FindNN(X, k * 10)
        for (i in 1:M) {
            j <- 1
            for (bla in 1:length(samp.points)) {
                n <- samp.points[bla]
                knnlen1 <- 0
                knnlen2 <- 0
                for (trial in 1:N) {
                  indices <- sample(dim(X)[1])
                  indices <- indices[1:n]
                  Dr <- D[[1]][indices, ]
                  Drr <- Dr[, indices]
                  L <- 0
                  Drr <- apply(Drr, 2, sort)
                  for (l in 1:dim(Drr)[2]) {
                    bla <- Drr[, l]
                    bla[bla == 0] <- NA
                    ind <- which.min(bla)
                    L = L + sum(Drr[(ind + 1):min(c(ind + k, 
                      dim(Drr)[2])), l])
                  }
                  knnlen1 <- knnlen1 + L
                  knnlen2 <- knnlen2 + L^2
                }
                knnlenavg.vec[i, j] <- knnlen1/N
                knnlenstd.vec[i, j] <- sqrt((knnlen2 - (knnlen1/N)^2 * 
                  N)/(N - 1))
                j <- j + 1
            }
            A <- cbind(log(samp.points), rep(1, Q))
            sol1 <- solve(t(A) \%*\% A) \%*\% t(A) \%*\% matrix(log(knnlenavg.vec[i, 
                ]), ncol = 1)
            dvec[i] <- gamma/(1 - sol1[1])
        }
        no.dims <- mean(abs(dvec))
    }
    if (method == "EigValue") {
        pca <- prcomp(X)
        mapped.x <- pca$x
        lambda <- pca$sdev^2/sum(pca$sdev^2)
        no.dims <- 0
        while (no.dims < dim(X)[2] - 1 & lambda[no.dims + 1] > 
            0.025) {
            no.dims <- no.dims + 1
        }
    }
    if (method == "MLE") {
        k1 <- 6
        k2 <- 12
        X <- t(X)
        d <- dim(X)[1]
        n <- dim(X)[2]
        X2 <- as.matrix(apply(X^2, 2, sum))
        knnmatrix <- matrix(0, k2, n)
        if (n < 3000) {
            distance <- t(repmat(X2, 1, n)) + repmat(X2, 1, n) - 
                2 * t(X) \%*\% X
            distance <- apply(distance, 1, sort)
            knnmatrix <- 0.5 * log(distance[2:(k2 + 1), ])
        }
        else {
            for (i in 1:n) {
                distance <- sort(repmat(as.matrix(X2[i]), 1, 
                  n) + t(X2) - 2 * t(X[, i]) \%*\% X)
                knnmatrix[, i] <- 0.5 * (log(distance[2:(k2 + 
                  1)]))
            }
        }
        S <- apply(knnmatrix, 2, cumsum)
        indexk <- repmat(as.matrix(k1:k2), 1, n)
        dhat <- -(indexk - 2)/(S[k1:k2, ] - knnmatrix[k1:k2, 
            ] * indexk)
        no.dims <- mean(dhat)
    }
    return(no.dims)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
