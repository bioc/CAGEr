#` plotAnnot
#'
#' Plot annotation statistics
#' 
#' Plot maping statistics of an object containing mapping statistics in counts as
#' percentages in stacked barplots.
#' 
#' Stacked barplots with error bars inspired from
#' <http://stackoverflow.com/questions/10417003/stacked-barplot-with-errorbars-using-ggplot2>.
#' See <http://www.biomedcentral.com/1471-2164/14/665/figure/F1> for example.
#' 
#' @param x An object from which can be extracted a table with columns named
#'   \code{promoter}, \code{exon}, \code{intron}, \code{mapped},
#'   \code{extracted}, \code{rdna}, and \code{tagdust}, that will be passed
#'   to the \code{mapStats} function.
#'   
#' @param scope The value on which to normalise (see the plotAnnot vignette).
#' 
#' @param title The title of the plot.
#' 
#' @param group A factor to group the samples.
#' 
#' @param customScope A function passed to the internal function \code{\link{mapStats}}
#'   for the definition of custom scopes.
#'   
#' @param normalise Whether to normalise or not. Default: TRUE.
#' 
#' @family CAGEr annotation functions
#' @family CAGEr plot functions
#' 
#' @author Charles Plessy
#' 
#' @examples
#' # ce <- To Do: prepare an example CAGEexp object
#' p <- plotAnnot(ce, 'qc', 'Here is the title')
#' print(p)
#' p + ggplot2::theme_bw()
#' ggplot2::theme_set(ggplot2::theme_bw()) ; p
#' plotAnnot(ce, 'qc', 'Same, non-normalised', normalise = FALSE)
#' 
#' @docType methods
#' @import ggplot2
#' @export

setGeneric("plotAnnot", function( x, scope, title, group = "default"
                                 , customScope = NULL, normalise = TRUE)
  standardGeneric("plotAnnot"))

setMethod("plotAnnot", "data.frame",
  function( x
          , scope
          , title
          , group = "default"
          , customScope = NULL
          , normalise = TRUE) {
  ggplot( mapStats(x, scope=scope, group=group, customScope = customScope, normalise = normalise)
        , aes( x    = group
             , y    = value
             , fill = variable)
        , main = all) +
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    geom_segment(aes( xend = group
                    , y    = ystart
                    , yend = yend)) +
    geom_point( aes( x = group
                   , y = yend)
              , shape = "|") +
    coord_flip() +
    ggtitle(title)
})

setMethod("plotAnnot", "DataFrame",
  function( x, scope, title, group, customScope, normalise) {
  plotAnnot( data.frame(x, check.names = FALSE)
           , scope=scope, title=title, group=group, customScope=customScope, normalise=normalise)
})

setMethod("plotAnnot", "CAGEexp",
  function( x, scope, title, group, customScope, normalise) {
  if (missing(group)) group <- sampleLabels(x)
  if (missing(title)) title <- paste("CAGEr object", dQuote(deparse(substitute(x))))
  plotAnnot( colData(x)
           , scope=scope, title=title, group=group, customScope=customScope, normalise=normalise)
})

#' mapStats
#' 
#' Process mapping statistics
#' 
#' Using a data frame containing mapping statistics in counts, transform the data in
#' percentages that can be used for stacked barplots.
#' 
#' See http://stackoverflow.com/questions/10417003/stacked-barplot-with-errorbars-using-ggplot2 about stacked barplot.
#' 
#' The \dQuote{mapped} and \dQuote{counts} scopes assume that transcript counts are available.
#' 
#' @param libs A data frame with columns named \code{promoter}, \code{exon}, \code{intron}
#'        \code{mapped}, \code{extracted}, \code{rdna}, and \code{tagdust}.
#' @param scope The value on which to normalise. \dQuote{all} and \dQuote{qc}
#'        normalise on the number of extracted tags, \dQuote{annotation} on the
#'        number of aligned tags, \dQuote{mapped} on the number of aligned tags
#'        and \dQuote{counts} on the transcript counts.
#' @param group A vector of factors defining groups in the data.  By default,
#'        the \dQuote{group} column of the \dQuote{libs} table.
#' @param customScope A function that implements a custom scope.  Use with
#'        \code{scope = "custom"}.  The function takes a data frame in input
#'        and returns a named list containing a data frame (\dQuote{libs}),
#'        a character vector of columns to be plotted (\dQuote{columns}), and
#'        a numeric vector of totals for the normalisation (\dQuote{total}).
#' @param normalise Whether to normalise or not. Default: TRUE.
#'
#' @return
#' Returns a data frame with mean and standard deviation of normalised mapping statistics,
#' plus absolute positions for the error bars.  The first column, \code{group}, is
#' a vector of factors sorted with the \code{gtools::mixedorder} function.
#' 
#' @details See the plotAnnot vignette for details on what the scopes are.
#' 
#' @author Charles Plessy
#' 
#' @family smallCAGEqc annotation functions
#'  
#' @examples
#' example(loadMoiraiStats)
#' mapStats(libs, "qc")
#' 
#' @importFrom gtools mixedorder
#' @importFrom reshape melt

mapStats <- function( libs
                    , scope = c( "all"
                               , "annotation"
                               , "counts"
                               , "mapped"
                               , "qc"
                               , "steps"
                               , "custom")
                    , group="default"
                    , customScope = NULL
                    , normalise = TRUE)
{
  scope <- match.arg(scope)
  if (identical(group, "default")) {
      if      ("group" %in% colnames(libs)) {
      group <- libs$group
    } else if ("Group" %in% colnames(libs)) {
      group <- libs$Group
    } else
      stop(paste("Missing", dQuote("group"), "column in the data frame."))
  }
  
  totalIs <- function(what) {
    if (normalise == FALSE) {
      total <<- 1
      return()
    }
    if (! what %in% colnames(libs))
      stop( paste(what, "column is missing, see the plotAnnot vignette.")
           , call. = FALSE)
    if (is.numeric(libs[, what])) {
      total <<- libs[, what]
    } else stop(paste(what, "column is not numeric."), call. = FALSE)
  }
    
  defaultToZero <- function(what)
    if (! (what %in% colnames(libs)))
      libs[, what] <<- 0
  
  defaultToZero("tagdust")
  
  if (scope == "counts") {
    totalIs("librarySizes")
    columns <- c("Promoter", "Exon", "Intron", "Intergenic")
    libs <- within(libs, {
      Promoter   <- promoter
      Exon       <- exon
      Intron     <- intron
      Intergenic <- librarySizes - promoter - intron - exon
    })
  } else if (scope == "custom") {
    stopifnot(is.function(customScope))
    custom.list <- customScope(libs)
    libs    <- custom.list$libs
    columns <- custom.list$columns
    total   <- custom.list$total
  } else if (scope == "mapped") {
    totalIs("mapped")
    columns <- c( "Promoter", "Exon", "Intron", "Intergenic"
                , "Duplicates", "Non_proper")
    libs <- within(libs, {
      Non_proper <- mapped - properpairs
      Duplicates <- properpairs - librarySizes
      Intergenic <- librarySizes - promoter - intron - exon
      Intron     <- intron
      Exon       <- exon
      Promoter   <- promoter
    })
  } else if (scope == "qc") {
    totalIs("extracted")
    columns <- c( "Tag_dust", "rDNA", "Spikes", "Unmapped"
                , "Non_proper", "Duplicates", "librarySizes")
    libs <- within(libs, {
      Tag_dust     <- extracted   - rdna - spikes - cleaned
      rDNA         <- rdna
      Spikes       <- spikes
      Unmapped     <- cleaned     - mapped
      Non_proper   <- mapped      - properpairs
      Duplicates   <- properpairs - librarySizes
      Counts       <- librarySizes
    })
   } else if (scope == "steps") {
    totalIs("extracted")
    columns <- c("Cleaning", "Mapping", "Deduplication", "librarySizes")
    libs <- within(libs, {
      Cleaning      <- extracted   - cleaned
      Mapping       <- cleaned     - properpairs
      Deduplication <- properpairs - librarySizes
      Counts        <- librarySizes
    })
    if ("total" %in% colnames(libs)) {
      totalIs("total")
      libs$Extraction <- with(libs, total - extracted)
      columns <- c("Extraction", columns)
    }
   } else {
    if (scope == "all")        totalIs("extracted")
    if (scope == "annotation") totalIs("mapped")
    columns <- c("promoter","exon","intron","mapped","rdna", "tagdust")
    libs <- within(libs, {
       mapped <- mapped - promoter - intron - exon
    })
  }
  
  doMean <- function (X) tapply(libs[,X] / total, group, mean)
  doSd   <- function (X) tapply(libs[,X] / total, group, sd  )
  
  # "simplify" needs to be FALSE so that conversion to data frame works even
  # when the group contains only a single level.
  mapstats          <- data.frame(sapply(columns, doMean, simplify = FALSE))
  mapstats$group    <- rownames(mapstats)
  mapstats[gtools::mixedorder(mapstats$group), ]
  mapstats$group    <- factor(mapstats$group, unique(mapstats$group))
  
  mapstats.sd       <- data.frame(sapply(columns, doSd, simplify = FALSE))
  mapstats.sd$group <- rownames(mapstats.sd)
  
  mapstats          <- reshape::melt(mapstats,    id.vars="group")
  mapstats$sd       <- reshape::melt(mapstats.sd, id.vars="group")$value
  
  mapstats          <- plyr::ddply( mapstats
                                  , plyr::.(group)
                                  , transform
                                  , ystart = cumsum(value)
                                  , yend   = cumsum(value) + sd)
  
  mapstats
}


#' @name annotateCTSS
#' 
#' @title Annotate and compute summary statistics.
#' 
#' @details Annotate the CTSS of a CAGEr object and compute annotation statistics.
#' 
#' @param object A \code{\link{CAGEexp}} object (\code{\link{CAGEset}}s are
#'   not supported).
#'   
#' @param ranges A \code{\link{GRanges}} object containing \code{gene_name},
#'   \code{type} and \code{transcript_type} metadata.
#'   
#' @return The input object with the following modifications:
#' 
#' \itemize{
#'   \item The Genomic Ranges of the \code{tagCountMatrix} experiment cains a
#'     \code{annotation} metadata column, with levels such as \dQuote{promoter},
#'     \dQuote{exon}, \dQuote{intron} and \dQuote{unknown"}.
#'   \item New column data added, indicating total counts in each of the annotation
#'     levels.
#' }
#' 
#' @seealso The \code{\link{Zv9_annot}} example data.
#' @family CAGEr object modifiers
#' @family CAGEr annotation functions
#' 
#' @author Charles Plessy
#' 
#' @examples 
#' \dontrun{
#' annotateCTSS(object)
#' }
#' 
#' @export

setGeneric("annotateCTSS", function(object, ranges) standardGeneric("annotateCTSS"))

setMethod("annotateCTSS", "CAGEset", function (object, ranges){
  stop("CAGEset objects not supported.")})

setMethod("annotateCTSS", "CAGEexp", function (object, ranges){
  objName <- deparse(substitute(object))
  if(is.null(experiments(object)$tagCountMatrix))
    stop(objName, " does not contain CTSS expressiond data, see ", dQuote("getCTSS()"), ".")
  
  CTSScoordinatesGR(object)$genes      <- ranges2genes(CTSScoordinatesGR(object), ranges)
  CTSScoordinatesGR(object)$annotation <- ranges2annot(CTSScoordinatesGR(object), ranges)
  
  annot <- sapply( CTSStagCountDF(object)
                 , function(X) tapply(X, CTSScoordinatesGR(object)$annotation, sum))
  colData(object)[levels(CTSScoordinatesGR(object)$annotation)] <- DataFrame(t(annot))
  
  if (validObject(object)) {
    assign(objName, object, envir = parent.frame())
    invisible(1)
  }
})

#' ranges2annot
#' 
#' hierarchical annotation of a Genomic Range as promoter, exon or intron.
#' 
#' @param ranges The Genomics Ranges object, for example extracted from a
#'               RangedSummarizedExperiment object with the \code{rowRanges}
#'               command.
#' 
#' @param annot A \code{\link{GRanges}} object containing \code{type} and
#'   \code{transcript_type} metadata providing enough information for
#'   inferring promoters, exons and introns.  Typically GENCODE.
#' 
#' @return A Run-length-encoded (Rle) factor of same length as the GRanges object,
#'         indicating if the interval is promoter, exon, intron or unknown.
#'         
#' @family CAGEr annotation functions
#' @seealso \code{\link{CTSScoordinatesGR}}, \code{\link{Zv9_annot}}
#' 
#' @author Charles Plessy
#' 
#' @examples
#' \dontrun{
#' # With a CAGEexp object called "ce" and a GRanges object called "gff"
#' # containing transcript/exon information like with GENCODE:
#' gff <- readRDS(system.file("extdata/Zv9_annot.rds", package = "CAGEr"))
#' CTSScoordinatesGR(ce)$annotation <- ranges2annot(CTSScoordinatesGR(ce), gff)
#' colData(ce)[levels(CTSScoordinatesGR(ce)$annotation)] <-
#'   DataFrame(t(sapply(CTSStagCountDF(ce), function(X) tapply(X, CTSScoordinatesGR(ce)$annotation, sum))))
#' }
#' 
#' @importFrom GenomicRanges findOverlaps promoters
#' @importFrom S4Vectors Rle

ranges2annot <- function(ranges, annot, showClasses=NULL) {
  if (is.null(annot$type) | is.null(annot$transcript_type))
    stop("Annotation must contain ", dQuote("type"), " and ", dQuote("transcript_type"), " metdata.")
  classes <- c("promoter", "exon", "intron", "unknown")
  typesWithPromoter <- c( "protein_coding", "processed_transcript", "lincRNA"
                          , "antisense", "processed_pseudogene"
                          , "unprocessed_pseudogene")
  
  if (! missing(showClasses) ) return(classes)
  
  findOverlapsBool <- function(A, B) {
    overlap <- findOverlaps(A, B)
    overlap <- as(overlap, "List")
    any(overlap)
  }
  
  p <- findOverlapsBool( ranges
                         , promoters( annot[ annot$type == "transcript"
                                           & annot$transcript_type %in% typesWithPromoter]
                                      , 500, 500))
  e <- findOverlapsBool(ranges, annot[annot$type == "exon"])
  t <- findOverlapsBool(ranges, annot[annot$type == "transcript"])
  
  annot <- sapply( 1:length(ranges), function(i) {
    if (p[i]) {classes[1]}
    else if (e[i]) {classes[2]}
    else if (t[i]) {classes[3]}
    else           {classes[4]}
  })
  annot <- factor(annot, levels = classes)
  Rle(annot)
}

#' ranges2genes
#' 
#' Assign gene symbol(s) to Genomic Ranges.
#' 
#' @param ranges Genomics Ranges object, for example extracted from a
#'               RangedSummarizedExperiment object with the \code{rowRanges}
#'               command.
#' 
#' @param genes A \code{\link{GRanges}} object containing \code{gene_name} metadata.
#' 
#' @return A character vector of same length as the GRanges object, indicating
#'         one gene symbol or a comma-separated list of gene symbols for each
#'         range.
#'         
#' @family CAGEr annotation functions
#' @seealso \code{\link{CTSScoordinatesGR}}, \code{\link{Zv9_annot}}
#' 
#' @author Charles Plessy
#' 
#' @examples
#' \dontrun{
#' # With a CAGEexp object called "ce" and a GRanges object called "gff"
#' # containing transcript/exon information like with GENCODE:
#' gff <- readRDS(system.file("extdata/Zv9_annot.rds", package = "CAGEr"))
#' CTSScoordinatesGR(ce)$genes <- ranges2genes(CTSScoordinatesGR(ce), gff)
#' }
#' 
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors List Rle unstrsplit
#' @importFrom IRanges extractList

ranges2genes <- function(ranges, genes) {
  if (is.null(genes$gene_name))
    stop("Annotation must contain ", dQuote("gene_name"), " metdata.")
  gnames <- findOverlaps(ranges, genes)
  gnames <- as(gnames, "List")
  gnames <- extractList(genes$gene_name, gnames)
  gnames <- unique(gnames)
  gnames <- unstrsplit(gnames, ";")
  Rle(gnames)
}

#' @name Zv9_annot
#' 
#' @title Example zebrafish annotation data
#' 
#' @details Annotation data for zebrafish's chromosome 17's interval  26000000-54000000
#' (Zv9/danRer7 genome), to be used in documentation examples.
#'
#' @docType data
#' @author Prepared by Charles Plessy \email{plessy@riken.jp} using archive ENSEMBL data.
#' @references \url{http://mar2015.archive.ensembl.org/biomart/}
#' @keywords data
#' 
#' @details Data was retreived from ENSEMBL's Biomart server using a query to extract
#' gene, transcripts and exon coordinates.  For the record, here it is as URL
#' (long, possibly overflowing).
#' 
#' http://mar2015.archive.ensembl.org/biomart/martview/78d86c1d6b4ef51568ba6d46f7d8b254?VIRTUALSCHEMANAME=default&ATTRIBUTES=drerio_gene_ensembl.default.structure.ensembl_gene_id|drerio_gene_ensembl.default.structure.ensembl_transcript_id|drerio_gene_ensembl.default.structure.start_position|drerio_gene_ensembl.default.structure.end_position|drerio_gene_ensembl.default.structure.transcript_start|drerio_gene_ensembl.default.structure.transcript_end|drerio_gene_ensembl.default.structure.strand|drerio_gene_ensembl.default.structure.chromosome_name|drerio_gene_ensembl.default.structure.external_gene_name|drerio_gene_ensembl.default.structure.gene_biotype|drerio_gene_ensembl.default.structure.exon_chrom_start|drerio_gene_ensembl.default.structure.exon_chrom_end|drerio_gene_ensembl.default.structure.is_constitutive|drerio_gene_ensembl.default.structure.rank&FILTERS=&VISIBLEPANEL=resultspanel
#' 
#' And here it is as XML.
#' 
#' \preformatted{<?xml version="1.0" encoding="UTF-8"?>
#' <!DOCTYPE Query>
#' <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
#'   <Dataset name = "drerio_gene_ensembl" interface = "default" >
#'     <Attribute name = "ensembl_gene_id" />
#'     <Attribute name = "ensembl_transcript_id" />
#'     <Attribute name = "start_position" />
#'     <Attribute name = "end_position" />
#'     <Attribute name = "transcript_start" />
#'     <Attribute name = "transcript_end" />
#'     <Attribute name = "strand" />
#'     <Attribute name = "chromosome_name" />
#'     <Attribute name = "external_gene_name" />
#'     <Attribute name = "gene_biotype" />
#'     <Attribute name = "exon_chrom_start" />
#'     <Attribute name = "exon_chrom_end" />
#'     <Attribute name = "is_constitutive" />
#'     <Attribute name = "rank" />
#'   </Dataset>
#' </Query>}
#' 
#' The downloaded file was then transformed as follows.
#' 
#' \preformatted{x <- read.delim("~/Downloads/mart_export.txt", stringsAsFactors = FALSE)
#' e <- GRanges(paste0("chr", x$Chromosome.Name), IRanges(x$Exon.Chr.Start..bp., x$Exon.Chr.End..bp.), ifelse(x$Strand + 1, "+", "-"))
#' e$gene_name <- Rle(x$Associated.Gene.Name)
#' e$transcript_type <- Rle(x$Gene.type)
#' e$type <- "exon"
#' e$type <- Rle(e$type)
#'
#' e <- GRanges(paste0("chr", x$Chromosome.Name), IRanges(x$Exon.Chr.Start..bp., x$Exon.Chr.End..bp.), ifelse(x$Strand + 1, "+", "-"))
#' e$gene_name <- Rle(x$Associated.Gene.Name)
#' e$transcript_type <- Rle(x$Gene.type)
#' e$type <- "exon"
#' e$type <- Rle(e$type)
#' e <- sort(unique(e))
#' 
#' g <- GRanges( paste0("chr", x$Chromosome.Name)
#'             , IRanges(x$Gene.Start..bp., x$Gene.End..bp.)
#'             , ifelse( x$Strand + 1, "+", "-"))
#'             
#' g$gene_name <- Rle(x$Associated.Gene.Name)
#' g$transcript_type <- Rle(x$Gene.type)
#' g$type <- "gene"
#' g$type <- Rle(g$type)
#' g <- sort(unique(g))
#' 
#' t <- GRanges( paste0("chr", x$Chromosome.Name)
#'             , IRanges(x$Transcript.Start..bp., x$Transcript.End..bp.)
#'             , ifelse( x$Strand + 1, "+", "-"))
#'             
#' t$gene_name <- Rle(x$Associated.Gene.Name)
#' t$transcript_type <- Rle(x$Gene.type)
#' t$type <- "transcript"
#' t$type <- Rle(t$type)
#' t <- sort(unique(t))
#' 
#' gff <- sort(c(g, t, e))
#' gff <- gff[seqnames(gff) == "chr17"]
#' gff <- gff[start(gff) > 26000000 & end(gff) < 54000000]
#' seqlevels(gff) <- seqlevelsInUse(gff)
#' 
#' saveRDS(gff, "inst/extdata/Zv9_annot.rds")}
#' 
#' @examples 
#' gff <- readRDS(system.file("extdata/Zv9_annot.rds", package = "CAGEr"))
NULL

#' @name CTSStoGenes
#' 
#' @title Make a gene expression table.
#' 
#' @description Add a gene expression table in the \code{GeneExpSE} experiment
#' of an annotated \code{\link{CAGEexp}} object.
#' 
#' @param object A CAGEexp object that was annotated with the \code{annotateCTSS} function.
#' 
#' @return The input object with the following modifications:
#' 
#' \itemize{
#'   \item A new \code{GeneExpSE} experiment containing gene expression levels as
#'     a \code{SummarizedExperiemnt} object with one assay called \code{counts}, which
#'     is a \code{DataFrame} of \code{Rle} integers.
#'   \item New \code{genes} column data added, indicating total number of gene symbols
#'     detected per library.
#' }
#' 
#' @author Charles Plessy
#' 
#' @seealso \code{\link{annotateCTSS}}.
#' 
#' @family CAGEr object modifiers
#' 
#' @examples 
#' \dontrun{
#' CTSStoGenes(ce)
#' }
#' 
#' @docType methods
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame Rle
#' @export

setGeneric("CTSStoGenes", function(object) standardGeneric("CTSStoGenes"))

setMethod( "CTSStoGenes"
         , signature(object = "CAGEset")
         , function (object)
  stop("Not supported for ", dQuote("CAGEset"), " objects.")
)

setMethod( "CTSStoGenes"
         , signature(object = "CAGEexp")
         , function (object){
  objName <- deparse(substitute(object))
  if (is.null(CTSScoordinatesGR(object)$genes))
    stop(objName, " is not annotated, see ", dQuote("annotateCTSS()"), ".")

  genes <- rowsum(CTSStagCountDf(object), as.factor(CTSScoordinatesGR(object)$genes))
  genes <- lapply(genes, Rle)
  genes <- DataFrame(genes, row.names = levels(factor(CTSScoordinatesGR(object)$genes)))

  GeneExpSE(object) <- SummarizedExperiment( assay   = SimpleList(counts = genes)
                                           , rowData = DataFrame(symbol = rownames(genes)))
  
  object$genes       <- colSums(as.data.frame(assay(GeneExpSE(object))) > 0)
  # object$geneSymbols <- countSymbols(assay(GeneExpSE(object)) %>% as.data.frame)
  
  if (validObject(object)) {
    assign(objName, object, envir = parent.frame())
    invisible(1)
  }
})