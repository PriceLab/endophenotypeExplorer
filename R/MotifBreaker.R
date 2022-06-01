library(R6)
library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
#--------------------------------------------------------------------------------
MotifBreaker = R6Class("MotifBreaker",
#--------------------------------------------------------------------------------
    private = list(
       targetGene = NULL,
       variants = NULL,
       motifs = NULL,
       quiet = NULL,
       results=NULL,
       scoreThreshold=NULL
       ),
#--------------------------------------------------------------------------------
    public = list(
       initialize = function(targetGene, variants, motifs, scoreThreshold, quiet=TRUE){
          private$targetGene <- targetGene
          private$variants <- variants
          private$motifs <- motifs
          private$quiet <- quiet
          }, # initialize

       prepareVariants = function(){
           if(!private$quiet)
               message(sprintf("starting lengthy snp lookup, creating GRanges from %d rsids",
                               length(private$variants)))
           elapsedTime <- 0
           if(class(private$variants) == "character"){  # need to create GRanges
             duration <- system.time(
                private$variants <- snps.from.rsid(rsid=private$variants,
                                                   dbSNP=SNPlocs.Hsapiens.dbSNP155.GRCh38,
                                                   search.genome=BSgenome.Hsapiens.UCSC.hg38))
             elapsedTime <- duration[["elapsed"]]
             }
           if(!private$quiet)
               message(sprintf("prepareVariants completing, GRanges of length %d in %5.2f seconds",
                               length(private$variants), elapsedTime))
          },

       getVariants = function(){
          private$variants
          },

       findBreaks = function(){
           worker.count <- min(length(private$variants), 40)
           if(!private$quiet)
               message(sprintf("worker.count: %d", worker.count))

           bpparam <- MulticoreParam(workers=worker.count)
           if(class(private$variants) == "character")
               self$prepareVariants()
           private$results <- motifbreakR(snpList = private$variants,
                                          filterp = TRUE,
                                          #threshold=0.0,
                                          pwmList = private$motifs,
                                          show.neutral=FALSE,
                                          method = c("ic", "log", "notrans")[1],
                                          bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                          BPPARAM = bpparam,
                                          verbose=!private$quiet)
           },

       getResults = function(){
           private$results
           }


      ) # public
#--------------------------------------------------------------------------------
) # MotifBreaker class

