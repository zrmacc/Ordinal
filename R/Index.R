#' Calculate the Probability Index.
#' 
#' Estimates the probability that the typical response in the 
#' treatment arm is greater than the typical response in the
#' reference arm, \eqn{P(Y_{1} > Y_{0})}.
#'
#' @param counts1 Counts per level in the tretment arm. 
#' @param counts0 Counts per level in the reference arm. 
#' @importFrom stats approxfun uniroot
#' @return Numeric interpolated quantile. 

ProbIndex <- function(counts1, counts0) {
  p1 <- counts1 / sum(counts1)
  p0 <- counts0 / sum(counts0)
  n <- length(p1)
  out <- 0
  for(i in 1:(n-1)){
    out <- out + sum(p0[i] * p1[(i+1):n])
  }
  return(out)
}


#' Generate Bootstrap Confidence Interval for the Difference in Probability Indices.
#'
#' Bootstrap confidence interval for the difference in probability indices:
#' \eqn{P(Y_{1} > Y_{0}) - P(Y_{0} > Y_{1})}.
#'
#' @param counts0 Counts per level in the reference arm. 
#' @param counts1 Counts per level in the tretment arm. 
#' @param B Bootstrap replicates.
#' @param alpha Alpha level.
#' @importFrom stats quantile 
#' @export
#' @return Numeric vector containing the probability indices for arms 1 `Arm1` and 0 `Arm0`, 
#' the difference `Delta`, and the lower `L` and upper `U` confidence bounds.

ProbIndexDiffCI <- function(counts0, counts1, B = 2000, alpha = 0.05) {
  
  # Observed indices.
  q0 <- ProbIndex(counts0, counts1)
  q1 <- ProbIndex(counts1, counts0)
  
  # Observed difference.
  dobs <- q1 - q0
  
  # Bootstrap function.
  aux <- function(b) {
    
    # Bootstrap data. 
    counts0b = BootCounts(counts0)
    counts1b = BootCounts(counts1)
    
    # Bootstrap quantiles.
    qb0 <- ProbIndex(counts0b, counts1b)
    qb1 <- ProbIndex(counts1b, counts0b)
    
    # Bootstrap difference.
    dboot <- qb1 - qb0
    return(dboot)
  }
  
  boot <- lapply(seq(1:B), aux)
  boot <- do.call(c, boot)
  
  # Confidence interval.
  alpha2 <- alpha / 2
  L <- as.numeric(quantile(boot, alpha2, na.rm = TRUE))
  U <- as.numeric(quantile(boot, 1 - alpha2, na.rm = TRUE))
  
  # Output.
  out <- c("Arm1" = q1, "Arm0" = q0, "Delta" = dobs, "L" = L, "U" = U)
  return(out)
}
  

#' Generate Permutation P-value for the Difference in in Probability Indices.
#'
#' Permutation p-value for the difference in probability indices:
#' \eqn{P(Y_{1} > Y_{0}) - P(Y_{0} > Y_{1})}.
#'
#' @param counts0 Counts per level in the reference arm. 
#' @param counts1 Counts per level in the tretment arm. 
#' @param q Quantile probability.
#' @param B Bootstrap replicates.
#' @importFrom stats quantile 
#' @export
#' @return Numeric vector containing the probability indices for arms 1 `Arm1` and 0 `Arm0`, 
#' the difference `Delta`, and the Pvalue. 

ProbIndexDiffP <- function(counts0, counts1, B = 2000, q = 0.5) {
  
  # Observed indices.
  q0 <- ProbIndex(counts0, counts1)
  q1 <- ProbIndex(counts1, counts0)
  
  # Observed difference.
  dobs <- q1 - q0
  
  # Function to bootstrap
  aux <- function(b) {
    
    # Permute treatment
    data_perm <- PermCounts(counts0, counts1)

    # Permuted quantiles. 
    qb0 <- ProbIndex(data_perm[, 2], data_perm[, 1])
    qb1 <- ProbIndex(data_perm[, 1], data_perm[, 2])
    
    # Bootstrap difference.
    dperm <- qb1 - qb0
    
    # Bootstrap difference as or more extreme
    out <- 1 * (abs(dperm) >= abs(dobs))
    return(out)
  }
  
  perm <- lapply(seq(1:B), aux)
  perm <- do.call(c, perm)
  
  # Permutation p-value.
  perm <- c(1, perm)
  p <- mean(perm)
  
  # Output.
  out <- c("Arm1" = q1, "Arm0" = q0, "Delta" = dobs, "P" = p)
  return(out)
}
  
  