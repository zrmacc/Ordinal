
# Description

This package provides functions for inference the difference in quantiles and probability indices comparing two treatment arms. 

# Installation


```r
devtools::install_github(repo = 'zrmacc/Ordinal')
```

# Examples

Simulated example data. 


```r
library(Ordinal)
arm0 <- as.numeric(rmultinom(n = 1, size = 200, prob = c(0.25, 0.25, 0.25, 0.25)))
arm1 <- as.numeric(rmultinom(n = 1, size = 200, prob = c(0.10, 0.30, 0.30, 0.30)))
names(arm0) = names(arm1) <- paste0('Category ', seq(1:4))
show(arm0)
```

```
## Category 1 Category 2 Category 3 Category 4 
##         43         45         62         50
```

```r
show(arm1)
```

```
## Category 1 Category 2 Category 3 Category 4 
##         28         46         70         56
```

## Difference in Quantiles

To find a confidence interval and p-vaue for the difference in medians $q = 0.5$:

```r
# Confidence interval calculation.
ci <- QuantDiffCI(
  counts0 = arm0,
  counts1 = arm1,
  q = 0.5
)
round(ci, digits = 2)
```

```
## Quantile     Arm1     Arm0    Delta        L        U 
##     0.50     2.82     2.33     0.49     0.25     0.71
```

```r
# P-value calculation. 
pval <- QuantDiffP(
  counts0 = arm0,
  counts1 = arm1,
  q = 0.5
)
round(pval, digits = 2)
```

```
## Quantile     Arm1     Arm0    Delta        P 
##     0.50     2.82     2.33     0.49     0.00
```

Note that the medians are linearly interpolated within integer bins. 

## Difference in Probability Indices

The difference in probability indices is defined as:
$$
\Delta = P(Y_{1} > Y_{0}) - P(Y_{0} > Y_{1}).
$$

To find a confidence interval and p-vaue for the difference in probability indices:

```r
# Confidence interval calculation.
ci <- ProbIndexDiffCI(
  counts0 = arm0,
  counts1 = arm1
)
round(ci, digits = 2)
```

```
##  Arm1  Arm0 Delta     L     U 
##  0.49  0.26  0.23  0.13  0.32
```

```r
# P-value calculation. 
pval <- ProbIndexDiffP(
  counts0 = arm0,
  counts1 = arm1
)
round(pval, digits = 2)
```

```
##  Arm1  Arm0 Delta     P 
##  0.49  0.26  0.23  0.00
```
