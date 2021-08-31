library(R6)
library(TrenaMultiScore)
library(TrenaProjectAD)
source("~/github/TrenaMultiScore/tools/runner/v2/tmsCore.R")
#--------------------------------------------------------------------------------
TrenaModelsSplitOnEQTL = R6Class("TrenaModelsSplitOnEQTL",
#--------------------------------------------------------------------------------
    private = list(
       targetGene = NULL,
       variant.rsid = NULL,
       tbl.fimo = NULL,
       tbl.oc   = NULL,
       etx = NULL,
       tmx = NULL,
       mtx = NULL
       ),
#--------------------------------------------------------------------------------
    public = list(
       initialize = function(targetGene, variant, tbl.fimo, tbl.oc, mtx){
          printf("initializing trenaModelsSplitOnEQTL")
          private$targetGene <- targetGene
          private$variant.rsid <- variant
          private$tbl.fimo <- tbl.fimo
          private$tbl.oc   <- tbl.oc
          private$mtx <- mtx
          private$etx <- EndophenotypeExplorer$new(targetGene, "hg38")
          private$tms <- TMS$new(tpad, targetGene, tbl.fimo=tbl.fimo, tbl.oc=tbl.oc)
          tms$addGeneHancer()
          tms$scoreFimoTFBS()
           },
       getTmsTable = function(){
          private$tms$getTfTable()
          }
       ) # public

) # class TrenaModelsSplitOnEQTL
#--------------------------------------------------------------------------------

