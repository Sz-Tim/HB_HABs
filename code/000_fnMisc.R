
# https://github.com/AndriSignorell/DescTools

ConDisPairs <-function(x){
  
  # tab is a matrix of counts
  # Based on code of Michael Friendly and Laura Thompson
  
  # slooooow because of 2 nested for clauses O(n^2)
  # this is NOT faster when implemented with a mapply(...)
  
  # Lookin for alternatives in C
  # http://en.verysource.com/code/1169955_1/kendl2.cpp.html
  # cor(..., "kendall") is for dimensions better
  
  n <- nrow(x)
  m <- ncol(x)
  pi.c <- pi.d <- matrix(0, nrow = n, ncol = m)
  
  row.x <- row(x)
  col.x <- col(x)
  
  for(i in 1:n){
    for(j in 1:m){
      pi.c[i, j] <- sum(x[row.x<i & col.x<j]) + sum(x[row.x>i & col.x>j])
      pi.d[i, j] <- sum(x[row.x<i & col.x>j]) + sum(x[row.x>i & col.x<j])
    }
  }
  C <- sum(pi.c * x)/2
  D <- sum(pi.d * x)/2
  
  return(list(pi.c = pi.c, pi.d = pi.d, C = C, D = D))
  
}

SomersDelta <- function(x,  y = NULL, direction=c("row","column"), conf.level = NA, ...) {
  
  if(!is.null(y)) tab <- table(x, y, ...)
  else tab <- as.table(x)
  
  # tab is a matrix of counts
  x <- ConDisPairs(tab)
  
  # use .DoCount
  #   if(is.na(conf.level)) {
  #     d.tab <- as.data.frame.table(tab)
  #     x <- .DoCount(d.tab[,1], d.tab[,2], d.tab[,3])
  #   } else {
  #     x <- ConDisPairs(tab)
  #   }
  
  m <- min(dim(tab))
  n <- sum(tab)
  switch( match.arg( arg = direction, choices = c("row","column") )
          , "row" = { ni. <- colSums(tab) }
          , "column" = { ni. <- rowSums(tab) }
  )
  wt <- n^2 - sum(ni.^2)
  # Asymptotic standard error: sqrt(sigma2)
  sigma2 <- 4/wt^4 * (sum(tab * (wt*(x$pi.c - x$pi.d) - 2*(x$C-x$D)*(n-ni.))^2))
  # debug: print(sqrt(sigma2))
  
  somers <- (x$C - x$D) / (n * (n-1) /2 - sum(ni. * (ni. - 1) /2 ))
  
  pr2 <- 1 - (1 - conf.level)/2
  ci <- qnorm(pr2) * sqrt(sigma2) * c(-1, 1) + somers
  
  if(is.na(conf.level)){
    result <- somers
  } else {
    result <- c(somers = somers,  lwr.ci=max(ci[1], -1), upr.ci=min(ci[2], 1))
  }
  
  return(result)
  
}
