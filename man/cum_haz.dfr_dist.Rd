% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dfr_dist.R
\name{cum_haz.dfr_dist}
\alias{cum_haz.dfr_dist}
\title{Method for obtaining the hazard function of a \code{dfr_dist} object.}
\usage{
\method{cum_haz}{dfr_dist}(x, ...)
}
\arguments{
\item{x}{The object to obtain the hazard function of.}

\item{...}{Additional arguments to pass into the \code{integrate} function.}
}
\value{
A function that computes the hazard function of the distribution.
It accepts \code{t}, the time at which to compute the hazard function, and
\code{par}, the parameters of the distribution. If \code{par} is \code{NULL}, then the
parameters of the \code{dfr_dist} object \code{x} are used. Finally, it passes any
additional arguments \code{...} to the \code{rate} function of the \code{dfr_dist}
object \code{x}.
}
\description{
Method for obtaining the hazard function of a \code{dfr_dist} object.
}
