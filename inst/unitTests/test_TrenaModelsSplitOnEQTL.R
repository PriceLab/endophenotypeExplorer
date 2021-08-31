library(EndophenotypeExplorer)
library(RUnit)
source("~/github/endophenotypeExplorer/R/TrenaModelsSplitOnEQTL.R")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))

   targetGene <- "PTK2B"
   variant <- "rs11780653"
   f <- "~/github/TrenaProjectAD/explore/ptk2b/tbl.fimo.PTK2B.RData"
   checkTrue(file.exists(f))
   tbl.fimo <- get(load(f))
   f <- "~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/mayoAllPeaks.1052789x4.RData"
   checkTrue(file.exists(f))
   tbl.atac <- get(load(f))
   splitter <- TrenaModelsSplitOnEQTL$new(targetGene, variant, tbl.fimo, tbl.atac)
   checkEquals(is(splitter), "TrenaModelsSplitOnEQTL")

} # test_ctor
#----------------------------------------------------------------------------------------------------
