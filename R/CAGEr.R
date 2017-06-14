#' @include AllClasses.R
#' @include CAGEexp.R
#' 
#' @name CAGEr-class
#' 
#' @title CAGEr objects
#' 
#' @details The CAGEr package provides two classes of objects to load, contain
#' and process CAGE data:
#' \itemize{
#'   \item The \code{\link{CAGEset}} class is the original object format in CAGEr,
#'   as when published in Haberle V, Forrest ARR, Hayashizaki Y, Carninci P and
#'   Lenhard B (2015). \dQuote{CAGEr: precise TSS data retrieval and high-resolution
#'   promoterome mining for integrative analyses.} \emph{Nucleic Acids Research},
#'   43, pp. e51., \href{http://nar.oxfordjournals.org/content/43/8/e51}
#'   {http://nar.oxfordjournals.org/content/43/8/e51}. 
#'   
#'   \item The \code{\link{CAGEexp}} class is a new class format in 2017, which
#'   is based on the \code{\link{MultiAssayExperiment}} class.  In comparison with
#'   \code{CAGEset}, objects, \code{CAGEexp} objects benefit from a a more efficient
#'   data storage, using \code{DataFrame}s of run-length-encodedy (\code{Rle})
#'   integers, allowing for the loading and use of much larger transcriptome datasets.
#' }
#' Most CAGEr functions support both classes interchangabely, and the \code{CAGEr}
#' class was created for methods that support both classes identically.
#' 
#' @docType class
#' @import methods
#' @exportClass CAGEr

setClassUnion("CAGEr", c("CAGEset", "CAGEexp"))