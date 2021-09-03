library(EndophenotypeExplorer)
library(TrenaProjectAD)
library(RUnit)
source("~/github/endophenotypeExplorer/R/TrenaModelsSplitOnEQTL.R")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_buildTable()
    test_buildModels()

} # runTests
#----------------------------------------------------------------------------------------------------
mayo.tcx.splitter <- function()
{
   targetGene <- "PTK2B"
   variant <- "rs11780653"

   f <- "~/github/TrenaProjectAD/explore/ptk2b/tbl.fimo.PTK2B.RData"
   checkTrue(file.exists(f))
   tbl.fimo <- get(load(f))

   f <- "~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/mayoAllPeaks.1052789x4.RData"
   checkTrue(file.exists(f))
   tbl.atac <- get(load(f))

   dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
   f <- "mtx.mayo.tcx.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x257.RData"
   full.path <- file.path(dir, f)
   checkTrue(file.exists(full.path))
   mtx.tcx <- get(load(full.path))
   checkEquals(dim(mtx.tcx), c(17009, 257))

   trenaProject <- TrenaProjectAD()

   splitter <- TrenaModelsSplitOnEQTL$new(targetGene, variant, tbl.fimo, tbl.atac,
                                          mtx.tcx, study.name="mayo", trenaProject)
   return(splitter)

} # mayo.tcx.splitter
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))

   splitter <- mayo.tcx.splitter()

   checkEquals(is(splitter), "TrenaModelsSplitOnEQTL")

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_buildTable <- function()
{
   message(sprintf("--- test_buildTable"))

   splitter <- mayo.tcx.splitter()
   splitter$buildGenomicTable()
   tbl.tms <- splitter$getTmsTable()
   checkEquals(dim(tbl.tms), c(51491, 17))

   tbl.sub <- subset(tbl.tms, abs(tss) < 1000 & oc & fimo_pvalue < 1e-7)
   tfs <- sort(unique(tbl.sub$tf))
   checkEquals(tfs, c("KLF6", "MAZ", "REST", "SP2", "ZNF263"))

} # test_buildTable
#----------------------------------------------------------------------------------------------------
test_buildModels <- function()
{
   message(sprintf("--- test_buildModels"))

   splitter <- mayo.tcx.splitter()
   splitter$buildGenomicTable()
   tbl.tms <- splitter$getTmsTable()
   checkEquals(dim(tbl.tms), c(51491, 17))

   tbl.filtered <- subset(tbl.tms, oc & fimo_pvalue < 1e-7)
   tbl.filtered  <- subset(tbl.tms, oc  & (chip | fimo_pvalue < 1e-4))
   tbl.filtered <- subset(tbl.tms, oc & grepl("PTK2B", annot.symbol))
   tfs <- sort(unique(tbl.filtered$tf))
   length(tfs)  # 513

   result <- splitter$buildModels(tbl.filtered)
   names(result)
   names(result$models)
   # use a value from colnames(result$models$trena.1)[-1]
   tbl.summary.rfScore <- splitter$summarizeStratifiedModels(result$models, "rfScore")
   tbl.summary.spearman <- splitter$summarizeStratifiedModels(result$models, "spearmanCoeff")
   tbl.summary.pearson <- splitter$summarizeStratifiedModels(result$models, "pearsonCoeff")
   tbl.summary.bicor <- splitter$summarizeStratifiedModels(result$models, "bicor")
   tbl.summary.xgboost <- splitter$summarizeStratifiedModels(result$models, "xgboost")
   tbl.summary.lasso <- splitter$summarizeStratifiedModels(result$models, "betaLasso")
   tbl.summary.ridge <- splitter$summarizeStratifiedModels(result$models, "betaRidge")

   tbl.all <- do.call(rbind,list(
                      tbl.summary.rfScore,
                      tbl.summary.spearman,
                      tbl.summary.pearson,
                      tbl.summary.bicor,
                      tbl.summary.xgboost,
                      tbl.summary.lasso,
                      tbl.summary.ridge))
   new.order <- order(tbl.all$tf)
   tbl.all <- tbl.all[new.order,]
   sort.order.tfs <- names(sort(table(tbl.all$tf), decreasing=TRUE))
   tbl.all <- tbl.all[subjectHits(findMatches(sort.order.tfs, tbl.all$tf)),]

} # test_buildModels
#----------------------------------------------------------------------------------------------------
viz <- function()
{
    if(!exists("igv"))
        igv <- start.igv("PTK2B", "hg38")
    zoomOut(); zoomOut()
    tbl.track <- data.frame(chrom="chr8", start=28046060-1, end=28046060, stringsAsFactors=FALSE)
    track <- DataFrameAnnotationTrack("rs11780653", tbl.track, color="red")
    displayTrack(igv, track)
    showGenomicRegion(igv, "chr8:27,230,585-28,124,907")

    tbl.track <- subset(tbl.atac, chrom=="chr8" & start >=27230585 & end <= 28124907)
    dim(tbl.track)
    track <- DataFrameQuantitativeTrack("atac", tbl.track, color="blue", autoscale=TRUE)
    displayTrack(igv, track)

    etx <- EndophenotypeExplorer$new("PTK2B", "hg38")
    tbl.eQTL <- etx$getEQTLsForGene()
    tbl.eQTL.sig <- subset(tbl.eQTL, pvalue < 0.01 & study=="ampad-mayo")
    dim(tbl.eQTL.sig)

    tbl.track <- tbl.eQTL.sig[, c("chrom", "hg38", "hg38", "pvalue")]
    colnames(tbl.track) <- c("chrom", "start", "end", "score")
    tbl.track$start <- tbl.track$start - 1
    tbl.track$score <- -log10(tbl.track$score)
    dim(tbl.track)
    track <- DataFrameQuantitativeTrack("eQTL", tbl.track, color="darkgreen", autoscale=TRUE)
    displayTrack(igv, track)

} # viz
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
