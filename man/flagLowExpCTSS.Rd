% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CAGEr.R
\name{flagLowExpCTSS}
\alias{flagLowExpCTSS}
\alias{flagLowExpCTSS,CAGEr-method}
\alias{flagLowExpCTSS,RangedSummarizedExperiment-method}
\alias{flagLowExpCTSS,DataFrame-method}
\alias{flagLowExpCTSS,matrix-method}
\alias{filterLowExpCTSS}
\alias{filterLowExpCTSS,CAGEr-method}
\title{Flag CTSSes based on sample expression}
\usage{
flagLowExpCTSS(
  object,
  threshold = 1,
  nrPassThreshold = 1,
  thresholdIsTpm = TRUE
)

\S4method{flagLowExpCTSS}{CAGEr}(
  object,
  threshold = 1,
  nrPassThreshold = 1,
  thresholdIsTpm = TRUE
)

\S4method{flagLowExpCTSS}{RangedSummarizedExperiment}(
  object,
  threshold = 1,
  nrPassThreshold = 1,
  thresholdIsTpm = TRUE
)

\S4method{flagLowExpCTSS}{DataFrame}(
  object,
  threshold = 1,
  nrPassThreshold = 1,
  thresholdIsTpm = TRUE
)

\S4method{flagLowExpCTSS}{matrix}(
  object,
  threshold = 1,
  nrPassThreshold = 1,
  thresholdIsTpm = TRUE
)

filterLowExpCTSS(
  object,
  threshold = 1,
  nrPassThreshold = 1,
  thresholdIsTpm = TRUE
)

\S4method{filterLowExpCTSS}{CAGEr}(
  object,
  threshold = 1,
  nrPassThreshold = 1,
  thresholdIsTpm = TRUE
)
}
\arguments{
\item{object}{An object from the \emph{CAGEr} package that contains expression
values for multiple samples.}

\item{threshold}{Flag CTSSs with signal \verb{< threshold}.}

\item{nrPassThreshold}{Only flag CTSSs when signal is below threshold in at
least \code{nrPassThreshold} samples.}

\item{thresholdIsTpm}{Logical, is threshold raw tag count value (\code{FALSE}) or
normalized signal (\code{TRUE}).}
}
\value{
\code{flagLowExpCTSS} returns a \code{\link{Rle}} vector where \code{TRUE} indicates the
index of a CTSS that passes the filter.

\code{filterLowExpCTSS} returns the \code{CAGEr} object where the output of
\code{flagLowExpCTSS} was stored internally.
}
\description{
Flag CTSSes for that do not pass an expression threshold in at least a given
number of samples.  This is typically used to ignore CTSSes that have been
seen only once in a single sample, as they can be considered to not be
reproduced.
}
\examples{
flagLowExpCTSS(exampleCAGEexp, threshold = 100, nrPassThreshold = 2)
}
\seealso{
Other CAGEr filter functions: 
\code{\link{filteredCTSSidx}()},
\code{\link{flagByUpstreamSequences}()}
}
\concept{CAGEr filter functions}
