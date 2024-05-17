#' @include CAGEexp.R CTSS.R
NULL

#' @noRd
#' 
#' @importFrom S4Vectors queryHits subjectHits runLength runValue
#' 
#' @examples 
#' # See also benchmarks/dominant_ctss.md
#' (ctss <- GRanges('chr1', IRanges(start = 1:10, end = 1:10), '+', score = c(1, 0, 0, 1, 2, 0, 2, 1, 0, 1)))
#' (clusters <- GRanges('chr1', IRanges(start = c(1,9), end = c(8,10)), '+'))
#' 
#' # The function assumes that all CTSSes have a score above zero
#' .ctss_summary_for_clusters(ctss[score(ctss)>0], clusters, removeSingletons = TRUE)
#' # If not the case, it will give incorrect nr_ctss and  fail to remove singletons
#' .ctss_summary_for_clusters(ctss, clusters, removeSingletons = TRUE)
#' 
#' # The function needs its output to be sorted and is not going to check it.
#' .ctss_summary_for_clusters(rev(ctss), clusters)
#' .ctss_summary_for_clusters(ctss, rev(clusters))
#' 
#' # Ties are resolved with 5' preference for both plus and minus strands.
#' # This may create a small bias.
#' .ctss_summary_for_clusters(ctss |> plyranges::mutate(strand = '-'), clusters |> plyranges::mutate(strand = '-'))

.ctss_summary_for_clusters <- function(ctss, clusters, removeSingletons = FALSE, keepSingletonsAbove = Inf) {
  # Match the clusters and the CTSS
  o <- findOverlaps(clusters, ctss)

  # number of CTSS per cluster
  rl <- rle(queryHits(o))$length
  
  # Where each run starts in the CTSS ranges object
  cluster_start_idx <- cumsum(c(1, head(rl, -1)))
  
  # Scores sorted by run and position
  grouped_scores <- extractList(score(ctss), o)
  
  # Find relative position of dominant CTSS in each run.
  # In case of ties, take the central one.  Note it might bias + and - strands
  # in a different way.
  local_max_idx <- sapply(grouped_scores, \(x) {
    w <- which(x == max(x))
    w[ceiling(length(w)/2)]
  })
  
  # Find absolute position of dominant CTSS in each run.
  global_max_ids <- cluster_start_idx + local_max_idx - 1
  
  # Record dominant CTSS as GRanges object.
  clusters$dominant_ctss <- granges(ctss)[subjectHits(o)][global_max_ids]
  
  # Record dominant CTSS score.  Mabye we should use its GRanges's score instead.
  clusters$tpm.dominant_ctss <-   score(ctss)[subjectHits(o)][global_max_ids]
  
  # Record total expression of the cluster
  score(clusters) <- Rle(sum(grouped_scores))

  # Count the number of clusters   
  clusters$nr_ctss <- rl
  
  # Remove clusters that match only one CTSS unless their expression is high enough
  if(removeSingletons)
    clusters <- subset(clusters, clusters$nr_ctss > 1 | score(clusters) >= keepSingletonsAbove)
  
  # Give numerical names to the clusters
  names(clusters) <- seq_along(clusters)
  
  # Finally, return the object
  clusters
}
