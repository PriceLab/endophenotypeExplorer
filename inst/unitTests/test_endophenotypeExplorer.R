library(RUnit)
library(EndophenotypeExplorer)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))
   vcf.base.url <- "https://igv-data.systemsbiology.net/static/ampad"
   vcf.data.file.chr2 <- "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_2.recalibrated_variants.vcf.gz"
   vcf.url <- sprintf("%s/%s", vcf.base.url, vcf.data.file.chr2)

   etx <- EndophenotypeExplorer$new("BIN1", "hg19", vcf.url)
   checkTrue(all(c("EndophenotypeExplorer", "R6") %in% class(etx)))

} # test_ctor
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()

