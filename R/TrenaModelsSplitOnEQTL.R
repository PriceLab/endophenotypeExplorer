library(R6)
library(TrenaMultiScore)
library(TrenaProjectAD)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
#--------------------------------------------------------------------------------
TrenaModelsSplitOnEQTL = R6Class("TrenaModelsSplitOnEQTL",
#--------------------------------------------------------------------------------
    private = list(
       study.name = NULL,
       targetGene = NULL,
       variant.rsid = NULL,
       tbl.fimo = NULL,
       tbl.oc   = NULL,
       etx = NULL,
       tms = NULL,
       mtx = NULL,
       tbl.tms.filtered = NULL   # filted by user, provides selected tfs
       ),
#--------------------------------------------------------------------------------
    public = list(
       initialize = function(targetGene, variant, tbl.fimo, tbl.oc, mtx,
                             study.name, trenaProject){
          stopifnot(study.name %in% c("rosmap", "sinai", "mayo"))
          private$study.name <- study.name
          private$targetGene <- targetGene
          private$variant.rsid <- variant
          private$tbl.fimo <- tbl.fimo
          private$tbl.oc   <- tbl.oc
          private$mtx <- mtx
          private$etx <- EndophenotypeExplorer$new(targetGene, "hg38")
          private$tms <- TMS$new(trenaProject, targetGene, tbl.fimo=tbl.fimo,
                                 tbl.oc=tbl.oc)
          },
       buildGenomicTable = function(){
          private$tms$addGeneHancer()
          private$tms$scoreFimoTFBS()
          },
       getTmsTable = function(){
          private$tms$getTfTable()
          },
       getEtx = function(){
          private$etx
          },
       buildModels = function(tbl.tms.filtered){
          tfs <- sort(unique(tbl.tms.filtered$tf))
          private$tbl.tms.filtered <- tbl.tms.filtered
          x <- private$etx$splitExpressionMatrixByMutationStatusAtRSID(private$mtx,
                                                                       private$variant.rsid,
                                                                       private$study.name)

          xx <-private$etx$trenaScoreGenotypeStratifiedExpression(x$wt, x$mut,
                                                                  private$targetGene, tfs)
          return(list(matrices=x, models=xx))
          }, # buildModels

        #------------------------------------------------------------
        # models is a list of two (wt/mut) trena models, and
        # a full set of wt vs mut tf deltas, one for each trena score:
        # lasso, spearmanCoeff, rfScore, etc.
        summarizeStratifiedModels = function(models, scoreName){

            stopifnot(scoreName %in% colnames(models$trena.1))
            tfs <- models[[scoreName]]$tf
            tbl.1 <- models$trena.1
            new.order <- order(abs(tbl.1[, scoreName]), decreasing=TRUE)
            tbl.1 <- tbl.1[new.order,]
            trena.wt.rank <- lapply(tfs,function(tf) grep(sprintf("^%s$", tf), tbl.1$gene))
            failures <- which(unlist(lapply(trena.wt.rank, length))==0)
            trena.wt.rank[failures] <- -1
            trena.wt.rank <- unlist(trena.wt.rank)

            tfs <- models[[scoreName]]$tf
            tbl.2 <- models$trena.2
            new.order <- order(abs(tbl.2[, scoreName]), decreasing=TRUE)
            tbl.2 <- tbl.2[new.order,]
            trena.mut.rank <- lapply(tfs, function(tf) grep(sprintf("^%s$", tf), tbl.2$gene))
            failures <- which(unlist(lapply(trena.mut.rank, length))==0)
            trena.mut.rank[failures] <- -1
            trena.mut.rank <- unlist(trena.mut.rank)

            models[[scoreName]]$wt.rank <- trena.wt.rank
            models[[scoreName]]$mut.rank <- trena.mut.rank
            threshold.delta <- fivenum(abs(models[[scoreName]]$delta))[4]
            tbl.delta.summary <- subset(models[[scoreName]],
                ((wt.rank > 0 & wt.rank <= 10) |
                 (mut.rank > 0 & mut.rank <= 10)) & abs(delta) > threshold.delta)
            if(nrow(tbl.delta.summary) > 0){
               tbl.delta.summary$rsid <- self$variant.rsid
               tf.counts <- unlist(lapply(tbl.delta.summary$tf,
                                          function(tf) length(grep(tf, private$tbl.tms.filtered$tf))))
               tbl.delta.summary$count <- tf.counts
               } # if nrow
           tbl.delta.summary
           } # summarizeStratifiedModels

       ) # public

) # class TrenaModelsSplitOnEQTL
#--------------------------------------------------------------------------------

