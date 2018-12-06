% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ars.R
\docType{class}
\name{namedVector}
\alias{namedVector}
\title{Class providing object with methods for communication with R6}
\format{\code{\link{R6Class}} object.}
\usage{
namedVector
}
\value{
Object of \code{\link{R6Class}} with methods to manipulate object
}
\description{
Class providing object with methods for communication with R6
}
\section{Fields}{

\describe{
\item{\code{names}}{Stores points}

\item{\code{values}}{Stores y values associated with points}
}}

\section{Methods}{

\describe{
  \item{\code{add(nams, vals)}}{This method adds names and values to the \code{namedVector}}

  \item{\code{delete_byname(nams)}}{This method deletes both the name and value of inputted nams from \code{namedVector}}
  \item{\code{sort_byname(decreasing=FALSE)}}{This method sorts the \code{namedVector} by name values.}}
}

\examples{
ex <- namedVector$new(names = c(1,2,3), values = c(4,5,6))
}
\keyword{datasets}