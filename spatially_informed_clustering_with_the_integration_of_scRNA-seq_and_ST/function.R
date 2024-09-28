# functions.R

library(matrixcalc)

# TopoLa function
TopoLa <- function(A, alpha) {
  A <- as.matrix(A)  
  A <- apply(A, 2, as.numeric)  
  temp <- A
  C <- solve((1/alpha * diag(ncol(temp)) + t(temp) %*% temp)) %*% t(temp) %*% temp
  Matrix <- A %*% C

  min_val <- min(Matrix)
  max_val <- max(Matrix)
  Matrix <- (Matrix - min_val) / (max_val - min_val)
  return(Matrix)
}

# Evaluation function
evalcluster <- function(truelabel, predlabel) {
  if(length(truelabel) != length(predlabel)) 
    stop("truelabel and predlabel must have the same length")
  
  total <- length(truelabel)
  x_ids <- unique(truelabel)
  y_ids <- unique(predlabel)
  
  # Mutual information
  MI <- 0.0
  for(idx in x_ids) {
    for(idy in y_ids) {
      idxOccur <- which(truelabel == idx)
      idyOccur <- which(predlabel == idy)
      idxyOccur <- intersect(idxOccur, idyOccur)
      if(length(idxyOccur) > 0) {
        MI <- MI + (length(idxyOccur) / total) * log2((length(idxyOccur) * total) / (length(idxOccur) * length(idyOccur)))
      }
    }
  }
  
  # Normalized Mutual information
  Hx <- 0
  for(idx in x_ids) {
    idxOccurCount <- length(which(truelabel == idx))
    Hx <- Hx - (idxOccurCount / total) * log2(idxOccurCount / total)
  }
  
  Hy <- 0
  for(idy in y_ids) {
    idyOccurCount <- length(which(predlabel == idy))
    Hy <- Hy - (idyOccurCount / total) * log2(idyOccurCount / total)
  }
  
  nmi <- 2 * MI / (Hx + Hy)
  
  # (adjusted) Rand Index
  tab <- table(truelabel, predlabel)
  conv_df <- as.data.frame.matrix(tab)
  n <- sum(tab)
  ni <- apply(tab, 1, sum)
  nj <- apply(tab, 2, sum)
  n2 <- choose(n, 2)
  nis2 <- sum(choose(ni[ni > 1], 2))
  njs2 <- sum(choose(nj[nj > 1], 2))
  ri <- 1 + (sum(tab^2) - (sum(ni^2) + sum(nj^2)) / 2) / n2
  ari <- (sum(choose(tab[tab > 1], 2)) - (nis2 * njs2) / n2) / ((nis2 + njs2) / 2 - (nis2 * njs2) / n2)
  
  out <- c(nmi, ri, ari)
  names(out) <- c("NMI", "RI", "ARI")
  return(out)
}

# K-means function
kmeansFunc_Iter <- function(data, k, centers_old = NULL) {
  set.seed(1234567)
  if(nrow(data) < 300000) {
    numStart <- 10
  } else {
    numStart <- 1
  }
  if(is.null(centers_old)) {
    cl <- suppressWarnings(try(kmeans(data, k, nstart = 1, iter.max = numStart), silent = TRUE))
  } else {
    cl <- suppressWarnings(try(kmeans(data, centers_old, nstart = 2, iter.max = numStart), silent = TRUE))
  }
  if(class(cl) == "try-error") {
    cl <- suppressWarnings(try(kmeans(data, k, nstart = 1, iter.max = numStart), silent = TRUE))
  }
  resList <- list(kmeans = cl$cluster, centers = cl$centers)
  return(resList)
}