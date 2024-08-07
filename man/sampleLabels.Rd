% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CAGEr.R, R/SetMethods.R
\name{sampleLabels}
\alias{sampleLabels}
\alias{sampleLabels,CAGEexp-method}
\alias{sampleLabels,CTSS-method}
\alias{sampleList}
\alias{sampleList,CAGEr-method}
\alias{sampleLabels<-}
\alias{sampleLabels<-,CAGEexp-method}
\alias{sampleLabels<-,CTSS-method}
\title{Get and set sample labels}
\usage{
sampleLabels(object)

\S4method{sampleLabels}{CAGEexp}(object)

\S4method{sampleLabels}{CTSS}(object)

sampleList(object)

\S4method{sampleList}{CAGEr}(object)

sampleLabels(object) <- value

\S4method{sampleLabels}{CAGEexp}(object) <- value

\S4method{sampleLabels}{CTSS}(object) <- value
}
\arguments{
\item{object}{A CAGEr object.}

\item{value}{A character vector with a unique and valid name for each sample.
The \code{names} attributes indicate the colors.}
}
\value{
\code{sampleLabels} returns a named character vector representing labels of all
CAGE datasets present in the CAGEr object.  The vector values are the labels and the
vector names are the colors.

\code{sampleList} returns a named list where elements and their
names are the sample names, for instance: \code{list(sampleA = "sampleA",
sampleB = "sampleB")}. Thus, after iterating on it with \code{lapply}, the
element names will be sample names.
}
\description{
\code{sampleLabels} gets or sets the labels and colors of CAGE datasets
(samples) from \code{\link{CAGEr}} objects.

\code{sampleList} is an accessory function for convenience
iteration in functions such as \code{\link{lapply}} or \code{\link{mapply}}.
There is no set method for \code{sampleList}.
}
\details{
In \code{CAGEexp} objects, renaming samples is possible only before
data is loaded.
}
\note{
If no colors are supplied, then default colors will be assigned
usign the \code{rainbow} function.  Assigned colors are not guaranteed
to be stable.
}
\examples{
sampleLabels(exampleCAGEexp)

sampleList(exampleCAGEexp)

}
\seealso{
\code{\link{setColors}}

Other CAGEr accessor methods: 
\code{\link{CTSScoordinatesGR}()},
\code{\link{CTSScumulativesTagClusters}()},
\code{\link{CTSSnormalizedTpmDF}()},
\code{\link{CTSStagCountDF}()},
\code{\link{GeneExpDESeq2}()},
\code{\link{GeneExpSE}()},
\code{\link{consensusClustersGR}()},
\code{\link{expressionClasses}()},
\code{\link{filteredCTSSidx}()},
\code{\link{genomeName}()},
\code{\link{inputFiles}()},
\code{\link{inputFilesType}()},
\code{\link{librarySizes}()},
\code{\link{seqNameTotalsSE}()},
\code{\link{tagClustersGR}()}

Other CAGEr setter methods: 
\code{\link{genomeName}()},
\code{\link{inputFiles}()},
\code{\link{inputFilesType}()},
\code{\link{setColors}()}
}
\author{
Vanja Haberle

Charles Plessy
}
\concept{CAGEr accessor methods}
\concept{CAGEr setter methods}
