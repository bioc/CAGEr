#' @include CAGEexp.R CTSS.R
NULL

#' Distance clustering
#' 
#' @param object The [`SummarizedExperiment::RangedSummarizedExperiment`] object
#'        containing CTSS information.
#' 
#' @param max.dist Maximal distance between two neighbouring CTSSs for them to
#'        be part of the same cluster.
#' 
#' @param removeSingletons Logical indicating if tag clusters containing only
#'        one CTSS be removed.
#' 
#' @param keepSingletonsAbove Controls which singleton tag clusters will be
#'        removed.  When `removeSingletons = TRUE`, only singletons with signal
#'        `< keepSingletonsAbove` will be removed.  Useful to prevent removing
#'        highly supported singleton tag clusters.  Default value `Inf` results
#'        in removing all singleton TCs when `removeSingletons = TRUE`.
#'        
#' @param removeSingletons Logical indicating if tag clusters containing only
#'        one CTSS be removed.
#' 
#' @param keepSingletonsAbove Controls which singleton tag clusters will be
#'        removed.  When `removeSingletons = TRUE`, only singletons with signal
#'        `< keepSingletonsAbove` will be removed.  Useful to prevent removing
#'        highly supported singleton tag clusters.  Default value `Inf` results
#'        in removing all singleton TCs when `removeSingletons = TRUE`.
#'        
#' @family CAGEr clustering methods
#' 
#' @examples 
#' distclu(CTSSnormalizedTpmGR(exampleCAGEexp, 1)[1:10])
#' distclu(CTSStagCountSE(exampleCAGEexp)[1:25,])
#' 
#' @export

setGeneric("distclu",
  function(
    object,
    max.dist = 20,
    removeSingletons = FALSE,
    keepSingletonsAbove = Inf
    ) standardGeneric("distclu"))

#' @rdname distclu

setMethod("distclu", "SummarizedExperiment",
          function(object, max.dist, removeSingletons, keepSingletonsAbove) {
  ctss.cluster.list <- GRangesList()
  for(s in colnames(object)) {
    message("\t-> ", s)
    d <- as(rowRanges(object), "CTSS")
    score(d) <- assays(object)[["normalizedTpmMatrix"]][[s]]
    d <- subset(d, score(d) > 0)
    ctss.cluster.list[[s]] <-
      distclu(d, max.dist = max.dist,
                removeSingletons = removeSingletons,
                keepSingletonsAbove = keepSingletonsAbove)
  }
  ctss.cluster.list
})

.distclu_CTSS <- function(object, max.dist, removeSingletons, keepSingletonsAbove) {
  clusters <- reduce(GRanges(object), min = max.dist)
  clusters <- .ctss_summary_for_clusters(object, clusters,
                                         removeSingletons    = removeSingletons,
                                         keepSingletonsAbove = keepSingletonsAbove)
  names(clusters) <- seq_along(clusters)
  clusters
}

#' @rdname distclu

setMethod("distclu", "CTSS", .distclu_CTSS)
