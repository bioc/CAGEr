#' @include CAGEexp.R CTSS.R
NULL

#' Distance clustering
#' 
#' The `"distclu"` method is an implementation of simple distance-based
#' clustering of data attached to sequences, where two neighbouring TSSs are
#' joined together if they are closer than some specified distance (see
#' [`GenomicRanges::reduce`] for implementation details.
#' 
#' Clustering is done for every CAGE dataset within the CAGEr object separately,
#' resulting in a different set of tag clusters for every CAGE dataset. TCs from
#' different datasets can further be aggregated into a single referent set of
#' consensus clusters by calling the [`aggregateTagClusters`] function.
#' 
#' @param object The [`SummarizedExperiment::RangedSummarizedExperiment`] object
#'        containing CTSS information, or just a [`CTSS`] object.
#' 
#' @param maxDist Maximal distance between two neighbouring CTSSs for them to
#'        be part of the same cluster.
#' 
#' @param keepSingletonsAbove Remove "singleton" tag clusters of width 1 with
#'        signal `< keepSingletonsAbove`.  Default value `0` results in keeping
#'        all TCs by default.  Setting it to `Inf` removes all singletons.
#'        
#' @family CAGEr clustering methods
#' @family CAGEr object modifiers
#' @family CAGEr clusters functions
#' 
#' @seealso [`aggregateTagClusters`]
#' 
#' @author Vanja Haberle
#' @author Charles Plessy
#' 
#' @return For `CTSS` input, a [`TagClusters`] object, for
#' `SummarizedExperiment` input, a [`GRangesList`] of [`TagClusters`] objects,
#' and for [`CAGEexp`] input, a modified object containing the tag clusters
#' stored as a `GRangesList` of [`TagClusters`] objects in its metadata slot
#' `tagClusters`.
#' 
#' @examples 
#' distclu(CTSSnormalizedTpmGR(exampleCAGEexp, 1)[1:10])
#' distclu(CTSStagCountSE(exampleCAGEexp)[1:25,])
#' ce <- distclu(exampleCAGEexp, maxDist = 20, keepSingletonsAbove = 100)
#' tagClustersGR(ce, "Zf.30p.dome")
#' 
#' @export

setGeneric("distclu",
  function(
    object,
    maxDist = 20,
    keepSingletonsAbove = 0
    ) standardGeneric("distclu"))

#' @rdname distclu

setMethod("distclu", "SummarizedExperiment",
          function(object, maxDist, keepSingletonsAbove) {
  ctss.cluster.list <- GRangesList()
  for(s in colnames(object)) {
    message("\t-> ", s)
    d <- as(rowRanges(object), "CTSS")
    score(d) <- assays(object)[["normalizedTpmMatrix"]][[s]]
    d <- subset(d, score(d) > 0)
    ctss.cluster.list[[s]] <-
      distclu(d, maxDist = maxDist,
                keepSingletonsAbove = keepSingletonsAbove)
  }
  endoapply(ctss.cluster.list, as, "TagClusters")
})

.distclu_CTSS <- function(object, maxDist, keepSingletonsAbove) {
  clusters <- reduce(GRanges(object), min.gapwidth = maxDist)
  clusters <- .ctss_summary_for_clusters(object, clusters)
  # Remove clusters that match only one CTSS unless their expression is high enough
  clusters <- subset(clusters, clusters$nr_ctss > 1 | score(clusters) >= keepSingletonsAbove)
  names(clusters) <- seq_along(clusters)
  as(clusters, "TagClusters")
}

#' @rdname distclu

setMethod("distclu", "CTSS", .distclu_CTSS)

#' @rdname distclu

setMethod( "distclu", "CAGEexp", function(object, maxDist, keepSingletonsAbove) {
  if (! "normalizedTpmMatrix" %in% assayNames(CTSStagCountSE(object)))
    stop( "Could not find normalized CAGE signal values, see ?normalizeTagCount.\n"
    , "distclu() needs normalized values to create its output tables, that "
    , "include TPM expression columns.")
  message("Clustering...")
  ctss.cluster.list <-
    distclu( object = CTSStagCountSE(object)[filteredCTSSidx(object),]
           , maxDist = maxDist, keepSingletonsAbove = keepSingletonsAbove)
  seqlevels(ctss.cluster.list) <- seqlevels(CTSStagCountSE(object))
  seqinfo(ctss.cluster.list)   <- seqinfo(CTSStagCountSE(object))
  # Changing the sequence levels may change the sort order.  Re-sort
  ctss.cluster.list <- sort(ctss.cluster.list)
  metadata(object)$tagClusters <- ctss.cluster.list
  object
})
