library(EndophenotypeExplorer)
library(TrenaProjectAD)
library(RUnit)
source("~/github/endophenotypeExplorer/R/MotifBreaker.R")
#----------------------------------------------------------------------------------------------------
fast.ndufs2.breaker <- function()
{
   targetGene <- "NDUFS2"
     #variants <- c("rs1136224","rs33941127","rs2070901","rs2070902","rs4233366")
     #variants <- get(load("../extdata/rsids-100-for-NDUFS2.RData"))
    variants <- get(load("gr.ndufs2.variants.RData"))
    tfs <- c("HSF2","RREB1","TGIF2","RXRG","NFIA","SOX9")
       #motifs <- query(MotifDb, c("jaspar2018", "sapiens")) # , tfs)
    motifs <- get(load("../extdata/motifs-7-for-NDUFS2.RData"))
    MotifBreaker$new(targetGene, variants, motifs, quiet=FALSE)

} # fast.ndufs2.breaker
#----------------------------------------------------------------------------------------------------
slow.ptk2b.breaker <- function()
{
   targetGene <- "PTK2B"
   motifs.selected <- query(MotifDb, c("hsapiens", "jaspar2018"))
   f <- "~/github/endophenotypeExplorer/inst/unitTests/motifBreakerData/eqtl.ptk2b.100.RData"
   get(load(f))  # two variables: 100 rsids in top.eqtl.rsids, GRanges length 128 in gr.variants
   MotifBreaker$new(targetGene,
                    variants=gr.variants, # top.eqtl.rsids,
                    motifs=motifs.selected,
                    quiet=FALSE)

} # slow.ptk2b.breaker
#----------------------------------------------------------------------------------------------------
if(!exists("breaker")){
    breaker <- fast.ndufs2.breaker()
   }
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    # test_prepareVariants()
    test_findBreaks.fast()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))
   checkEquals(is(breaker), "MotifBreaker")

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_prepareVariants <- function()
{
   message(sprintf("--- test_prepareVariants"))

   targetGene <- "NDUFS2"
   rsids <- c("rs1136224","rs33941127","rs2070901","rs2070902","rs4233366")
   tfs <- c("HSF2","RREB1","TGIF2","RXRG","NFIA","SOX9")
     #motifs <- query(MotifDb, c("jaspar2018", "sapiens")) # , tfs)
   motifs <- get(load("../extdata/motifs-7-for-NDUFS2.RData"))
   breaker <- MotifBreaker$new(targetGene, rsids, motifs, quiet=FALSE)

   char.variants <- breaker$getVariants()
   checkEquals(char.variants, rsids)

   breaker$prepareVariants()
   gr.variants <- breaker$getVariants()
   checkEquals(class(gr.variants), "GRanges")
   checkEquals(sort(names(gr.variants)), c("rs1136224","rs2070901","rs2070902","rs33941127","rs4233366"))

   checkEquals(as.character(gr.variants$REF), c("A","G","C","C","C"))
   checkEquals(as.character(gr.variants$ALT), c("G","T","T","T","T"))

} # test_findBreaks
#----------------------------------------------------------------------------------------------------
test_findBreaks.fast <- function()
{
    message(sprintf("--- test_findBreaks"))

    results <- breaker$findBreaks()
    tbl.results <- as.data.frame(results, row.names=NULL)
    lapply(tbl.results, class)
    colnames(tbl.results)[1] <- "chrom"
    tbl.results$chrom <- as.character(tbl.results$chrom)
    tbl.results <- subset(tbl.results, effect=="strong")
    checkEquals(nrow(tbl.results), 16)

} # test_findBreaks.fast
#----------------------------------------------------------------------------------------------------
test_runFull.ptk2b <- function()
{
    message(sprintf("--- test_runFull.ptk2b"))
    breaker <- slow.ptk2b.breaker()
    results <- breaker$findBreaks()
    save(results, file="ptk2b.top100eqtls.allMotifs.RData")
    

} # test_runFull.ptk2b
#----------------------------------------------------------------------------------------------------
# use NDUFS2 as a test case.  this scratchpad function is where eQTLs and candidate TFs were
# figured out
prep_NDUFS2 <- function()
{
   targetGene <- "NDUFS2"
   etx <- EndophenotypeExplorer$new(targetGene, "hg38")
   tbl.eQTL <- etx$getEQTLsForGene()
   tbl.eQTL <- tbl.eQTL[order(tbl.eQTL$pvalue, decreasing=FALSE),]
   dim(tbl.eQTL)
   library(trena)
   library(MotifDb)
   motifs <- query(MotifDb, c("sapiens", "jaspar2018"))
   f <- "~/github/TrenaProjectAD/inst/extdata/expression/temporalCortex.15167x264.RData"
   file.exists(f)
   mtx.rna <- get(load(f))
   tfs <- intersect(mcols(motifs)$geneSymbol, rownames(mtx.rna))
   length(tfs)

   solver <- EnsembleSolver(mtx.rna,
                            targetGene=targetGene,
                            candidateRegulators=tfs,
                            solverNames=c("lasso", "Ridge", "Spearman", "Pearson",
                                          "RandomForest", "xgboost", "bicor"))
   tbl.out <- run(solver)
   tbl.out <- tbl.out[order(abs(tbl.out$bicor), decreasing=TRUE),]
   tfs.oi <- head(tbl.out$gene)
   motifs <- query(motifs, andStrings="sapiens", orString=tfs.oi)
   save(motifs, file="../extdata/motifs-7-for-NDUFS2.RData")
   rsids <- head(tbl.eQTL$rsid, n=100)
   length(rsids)
   save(rsids, file="../extdata/rsids-100-for-NDUFS2.RData")

} # prep_NDUFS
#----------------------------------------------------------------------------------------------------
if(!interactive())
    test_runFull.ptk2b()

