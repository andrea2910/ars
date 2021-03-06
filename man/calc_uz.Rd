% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calc_uz}
\alias{calc_uz}
\title{calc_uz}
\usage{
calc_uz(y_l, y_r, hy_l, hy_r, hyprim_l, hyprim_r)
}
\arguments{
\item{y_l}{y points to the left}

\item{y_r}{y points to the right}

\item{hy_l}{h(y_l) values}

\item{hy_r}{h(y_r) values}

\item{hyprim_l}{h'(y_l)}

\item{hyprim_r}{h'(y_r)}
}
\value{
a namedVector object with names=z and values=uz
}
\description{
calculate the z points (or intersection between the tangent lines at given points) given a (sorted) set of points and their h values and h_prim values
}
