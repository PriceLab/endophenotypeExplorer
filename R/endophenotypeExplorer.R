# EndophenotypeExplorer
#----------------------------------------------------------------------------------------------------
#' @importFrom RPostgreSQL dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @import RPostgreSQL
#' @import org.Hs.eg.db
#' @import R6
#' @import VariantAnnotation
#' @import GenomicRanges
#' @importFrom rtracklayer liftOver import.chain
#'
#' @title EndophenotypeExplorer
#------------------------------------------------------------------------------------------------------------------------
#' @name EndophenotypeExplorer-class
#' @rdname EndophenotypeExplorer-class
#' @aliases EndophenotypeExplorer
#'
# library(R6)
# library(VariantAnnotation)
# library(RPostgreSQL)
#----------------------------------------------------------------------------------------------------
#' R6 Class for exploring associations between eQTLs, variants, gene expression and gene regulation
#'
#' An endophenotypeExplore has access to eQTLs, variant calls, and more
#'
#' @export
EndophenotypeExplorer = R6Class("EndophenotypeExplorer",
    #--------------------------------------------------------------------------------
    private = list(target.gene=NULL,
                   vcf.url=NULL,
                   geneReg.db=NULL,
                   default.genome=NULL
                   ),
    #--------------------------------------------------------------------------------
    public = list(
      #' @description
      #' Create a new EndophenotypeExplorer object
      #' @param target.gene  Gene of interest.
      #' @param default.genome UCSC code, either `hg19` or `hg38`.
      #' @param vcf.url https endpoint from serving indexed vcf files
      #' @return A new `EndophenotypeExplorer` object.

        initialize = function(target.gene, default.genome, vcf.url){
            private$target.gene <- target.gene
            private$default.genome <- default.genome
            private$vcf.url <- vcf.url
           }
       ) # public

    ) # class EndophenotypeExplorer

#--------------------------------------------------------------------------------


