#' Bootstrap counts. 
#'
#' @param counts Counts in each level.
#' @importFrom stats rmultinom
#' @return Bootstrapped counts. 

BootCounts <- function(counts) {
  n <- sum(counts)
  frq <- counts / n
  out <- as.numeric(rmultinom(n = 1, size = n, prob = frq))
  return(out)
}


#' Permute counts. 
#'
#' @param counts0 Counts per level in the reference arm. 
#' @param counts1 Counts per level in the tretment arm. 
#' @importFrom stats rmultinom
#' @return Levels by 2 matrix of permuted counts. 

PermCounts <- function(counts0, counts1) {
  n0 <- sum(counts0)
  n1 <- sum(counts1)
  frq <- (counts0 + counts1) / (n0 + n1)
  out0 <- as.numeric(rmultinom(n = 1, size = n0, prob = frq))
  out1 <- as.numeric(rmultinom(n = 1, size = n1, prob = frq))
  out <- cbind(out0, out1)
  colnames(out) <- c('Arm0', 'Arm1')
  return(out)
}
