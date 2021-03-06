% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Index.R
\name{ProbIndexDiffCI}
\alias{ProbIndexDiffCI}
\title{Generate Bootstrap Confidence Interval for the Difference in Probability Indices.}
\usage{
ProbIndexDiffCI(counts0, counts1, B = 2000, alpha = 0.05)
}
\arguments{
\item{counts0}{Counts per level in the reference arm.}

\item{counts1}{Counts per level in the tretment arm.}

\item{B}{Bootstrap replicates.}

\item{alpha}{Alpha level.}
}
\value{
Numeric vector containing the probability indices for arms 1 `Arm1` and 0 `Arm0`, 
the difference `Delta`, and the lower `L` and upper `U` confidence bounds.
}
\description{
Bootstrap confidence interval for the difference in probability indices:
\eqn{P(Y_{1} > Y_{0}) - P(Y_{0} > Y_{1})}.
}
