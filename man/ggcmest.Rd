% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cmest.R
\name{ggcmest}
\alias{ggcmest}
\title{Plotting Point Estimates and Confidence Intervals of Causal Effects}
\usage{
ggcmest(
  x,
  errorbar.width = 0.3,
  errorbar.size = 0.3,
  errorbar.colour = "black",
  point.size = 1,
  point.colour = "blue",
  refline = TRUE,
  refline.colour = "red",
  refline.size = 0.3
)
}
\arguments{
\item{x}{an object of class \code{cmest}}

\item{errorbar.width}{width of errorbars for confidence intervals. Default is \code{0.3}.}

\item{errorbar.size}{size of errorbars for confidence intervals. Default is \code{0.3}.}

\item{errorbar.colour}{colour of errorbars for confidence intervals. Default is \code{black}.}

\item{point.size}{size of points for point estimates. Default is \code{1}.}

\item{point.colour}{colour of points for point estimates. Default is \code{blue}.}

\item{refline}{a logical value. If \code{true}, include a reference line at
\code{y = 0} when effects are on the difference scale and include a reference line at
\code{y = 1} when effects are on the ratio scale. Default is \code{TRUE}.}

\item{refline.colour}{colour of the reference line. Default is \code{red}.}

\item{refline.size}{size of the reference line. Default is \code{0.3}.}
}
\description{
\code{ggcmest} is used to plot results of \code{cmest} nicely with plotting functions
in the \link{ggplot2} package. Additional layers can be added to this plot using other
plotting functions in the \link{ggplot2} package.
}
\examples{

library(CMAverse)
library(ggplot2)

x <- cmest(data = cma2020, model = "rb", outcome = "contY",
exposure = "A", mediator = "M2", basec = c("C1", "C2"),
EMint = TRUE, mreg = list("multinomial"), yreg = "linear",
astar = 0, a = 1, mval = list("M2_0"), estimation = "paramfunc",
inference = "delta")

ggcmest(x) +
theme(axis.text.x = element_text(angle = 45))

ggcmest(x) +
coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

}
\seealso{
\code{\link{cmest}}, \code{\link{ggplot2}}.
}
