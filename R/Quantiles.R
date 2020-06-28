#' Find Quantile of Ordinal Categorical Data. 
#'
#' @param counts Counts in each level. 
#' @param q Quantile.
#' @importFrom stats approxfun uniroot
#' @return Numeric interpolated quantile. 

FindQuantile <- function(counts, q) {
  offset_levels <- seq(from = 0.5, to = length(counts) + 0.5, by = 1)
  frq <- counts / sum(counts)
  offset_probs <- c(0, round(cumsum(frq), digits = 16))
  g <- approxfun(x = offset_levels, y = offset_probs)
  h <- function(x){g(x) - q}
  zero <- uniroot(
    f = h, 
    lower = min(offset_levels),
    upper = max(offset_levels)
  )
  return(zero$root)
}


#' Generate Bootstrap Confidence Interval for the Difference in Quantiles.
#'
#' Bootstrap confidence interval for the difference in quantiles.
#'
#' @param counts0 Counts per level in the reference arm. 
#' @param counts1 Counts per level in the tretment arm. 
#' @param q Quantile probability.
#' @param B Bootstrap replicates.
#' @param alpha Alpha level.
#' @importFrom stats quantile 
#' @export
#' @return Numeric vector containing the quantile probability `Quantile`, the
#'   estimated quantiles in arms 1 `Arm1` and 0 `Arm0`, the difference `Delta`,
#'   and the lower `L` and upper `U` confidence bounds.

QuantDiffCI <- function(counts0, counts1, q = 0.5, B = 2000, alpha = 0.05) {
  
  # Observed quantiles.
  q0 <- FindQuantile(counts0, q = q)
  q1 <- FindQuantile(counts1, q = q)
  
  # Observed difference.
  dobs <- q1 - q0
  
  # Bootstrap function.
  aux <- function(b) {
    
    # Bootstrap quantiles.
    qb0 <- FindQuantile(BootCounts(counts0), q = q)
    qb1 <- FindQuantile(BootCounts(counts1), q = q)
    
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
  out <- c("Quantile" = q, "Arm1" = q1, "Arm0" = q0, "Delta" = dobs, "L" = L, "U" = U)
  return(out)
}
  

#' Generate Permutation P-value for the Difference in Quantiles.
#'
#' Permutation p-value for the difference in quantiles.
#'
#' @param counts0 Counts per level in the reference arm. 
#' @param counts1 Counts per level in the tretment arm. 
#' @param q Quantile probability.
#' @param B Bootstrap replicates.
#' @importFrom stats quantile 
#' @export
#' @return Numeric vector containing the quantile probability `Quantile`, the
#'   estimated quantiles in arms 1 `Arm1` and 0 `Arm0`, the difference `Delta`,
#'   and the permutation p-value `P`.

QuantDiffP <- function(counts0, counts1, B = 2000, q = 0.5) {
  
  # Observed quantiles.
  q0 <- FindQuantile(counts0, q = q)
  q1 <- FindQuantile(counts1, q = q)
  
  # Observed difference.
  dobs <- q1 - q0
  
  # Function to bootstrap
  aux <- function(b) {
    
    # Permute treatment
    data_perm <- PermCounts(counts0, counts1)

    # Permuted quantiles. 
    qb0 <- FindQuantile(data_perm[, 1], q = q)
    qb1 <- FindQuantile(data_perm[, 2], q = q)
    
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
  out <- c("Quantile" = q, "Arm1" = q1, "Arm0" = q0, "Delta" = dobs, "P" = p)
  return(out)
}
  
  
  