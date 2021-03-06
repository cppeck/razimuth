% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{update_tuning}
\alias{update_tuning}
\title{Update tuning parameters}
\usage{
update_tuning(k, accept_tmp, tune)
}
\arguments{
\item{k}{iteration}

\item{accept_tmp}{acceptance rate during past 50 iterations}

\item{tune}{previous batch tuning parameter values}
}
\description{
Tuning parameter update during burn-in.
}
