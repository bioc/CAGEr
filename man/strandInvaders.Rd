% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StrandInvaders.R
\name{Strand invaders}
\alias{Strand invaders}
\alias{findStrandInvaders}
\alias{removeStrandInvaders}
\alias{findStrandInvaders,CAGEexp-method}
\alias{removeStrandInvaders,CAGEexp-method}
\alias{findStrandInvaders,CTSS-method}
\alias{removeStrandInvaders,CTSS-method}
\title{Detect and remove strand invasion artefacts}
\usage{
findStrandInvaders(object, distance = 1, barcode = NULL, linker = "TATAGGG")

removeStrandInvaders(object, distance = 1, barcode = NULL, linker = "TATAGGG")

\S4method{findStrandInvaders}{CAGEexp}(object, distance = 1, barcode = NULL, linker = "TATAGGG")

\S4method{removeStrandInvaders}{CAGEexp}(object, distance = 1, barcode = NULL, linker = "TATAGGG")

\S4method{findStrandInvaders}{CTSS}(object, distance = 1, barcode = NULL, linker = "TATAGGG")

\S4method{removeStrandInvaders}{CTSS}(object, distance = 1, barcode = NULL, linker = "TATAGGG")
}
\arguments{
\item{object}{A \code{\link{CAGEexp}} object object containing CTSS data and the name
of a reference genome.}

\item{distance}{The maximal edit distance between the genome and linker
sequences.  Regardless this parameter, only a single mismatch is
allowed in the last three bases of the linker.}

\item{barcode}{A vector of sample barcode sequences, or the name of a column
metadata of the \code{CAGEexp} object containing this information.
(\emph{Not implemented yet})}

\item{linker}{The sequence of the tail of the template-switching
oligonucleotide, that will be matched with the genome sequence
(defaults to \code{TATAGGG}).}
}
\value{
\code{findStrandInvaders} returns a logical-\link{Rle} vector indicating the
position of the strand invaders in the input ranges.

With \link{CTSS} objects as input \code{removeStrandInvaders} returns the
object after removing the CTSS positions identified as strand invaders.
In the case of \code{CAGEexp} objects, a modified object is returned.  Its sample
metadata is also updated by creating a new \code{strandInvaders} column that
indicates the number of molecule counts removed.  This value is subtracted
from the \code{counts} colum so that the total number of tags is still equal to
\code{librarySizes}.
}
\description{
\code{findStrandInvaders} detects strand invasion artefacts in the
CTSS data.  \code{removeStrandInvaders} removes them.

\emph{Strand invaders} are artefacts produced by \emph{template switching} reactions
used in methods such as \emph{nanoCAGE} and its derivatives (\emph{C1 CAGE}, ...).
They are described in details in Tang \emph{et al.}, 2013.  Briefly, these
artefacts create CAGE-like signal downstream of genome sequences highly
similar to the tail of template-switching oligonucleotides, which is
\code{TATAGGG} in recent (2017) nanoCAGE protocols.  Since these artefacts
represent truncated cDNAs, they do not indicate promoter regions.  It is
therefore advisable to remove these artefacts.  Moreover, when a sample
barcode is near the linker sequence (which is not the case in recent
nanoCAGE protocols), the strand-invasion artefacts can produce
\emph{sample-specific biases}, which can be confounded with biological effects
depending on how the barcode sequences were chosen.  A \code{barcode} parameter
is provided to incorporate this information.
}
\examples{
# Note that these examples do not do much on the example data since it was
# not constructed using a protocol based using the template-switching method.

findStrandInvaders(exampleCAGEexp)
removeStrandInvaders(exampleCAGEexp)

}
\references{
Tang \emph{et al.}, \dQuote{Suppression of artifacts and barcode bias in
high-throughput transcriptome analyses utilizing template switching.}
\emph{Nucleic Acids Res.} \strong{2013} Feb 1;41(3):e44.
PubMed ID: \href{https://pubmed.gov/23180801}{23180801},
DOI: \href{https://doi.org/10.1093/nar/gks1128}{10.1093/nar/gks112}
}
