#' @include CAGEexp.R CTSS.R
NULL

#' Distance clustering
#' 
#' @param object The [`SummarizedExperiment::RangedSummarizedExperiment`] object
#'        containing CTSS information, or just a [`CTSS`] object.
#' 
#' @param max.dist Maximal distance between two neighbouring CTSSs for them to
#'        be part of the same cluster.
#' 
#' @param keepSingletonsAbove Remove "singleton" tag clusters of width 1 with
#'        signal `< keepSingletonsAbove`.  Default value `0` results in keeping
#'        all TCs by default.  Setting it to `Inf` removes all singletons.
#'        
#' @family CAGEr clustering methods
#' 
#' @return For `CTSS` input, a [`TagClusters`] object, and for
#' `SummarizedExperiment` in put, a [`GRangesList`] of [`TagClusters`] objects.
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
    keepSingletonsAbove = 0
    ) standardGeneric("distclu"))

#' @rdname distclu

setMethod("distclu", "SummarizedExperiment",
          function(object, max.dist, keepSingletonsAbove) {
  ctss.cluster.list <- GRangesList()
  for(s in colnames(object)) {
    message("\t-> ", s)
    d <- as(rowRanges(object), "CTSS")
    score(d) <- assays(object)[["normalizedTpmMatrix"]][[s]]
    d <- subset(d, score(d) > 0)
    ctss.cluster.list[[s]] <-
      distclu(d, max.dist = max.dist,
                keepSingletonsAbove = keepSingletonsAbove)
  }
  endoapply(ctss.cluster.list, as, "TagClusters")
})

.distclu_CTSS <- function(object, max.dist, keepSingletonsAbove) {
  clusters <- reduce(GRanges(object), min = max.dist)
  clusters <- .ctss_summary_for_clusters(object, clusters,
                                         keepSingletonsAbove = keepSingletonsAbove)
  names(clusters) <- seq_along(clusters)
  as(clusters, "TagClusters")
}

#' @rdname distclu

setMethod("distclu", "CTSS", .distclu_CTSS)
