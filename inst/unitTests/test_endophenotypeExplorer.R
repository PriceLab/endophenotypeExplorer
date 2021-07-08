library(RUnit)
library(EndophenotypeExplorer)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_vcf.sampleID.to.clinicalTable()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   checkTrue(all(c("EndophenotypeExplorer", "R6") %in% class(etx)))
   expected <- "https://igv-data.systemsbiology.net/static/ampad/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_recalibrated_variants/chr2.vcf.gz"
   checkEquals(etx$getVcfUrl(), expected)

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_vcf.sampleID.to.clinicalTable <- function()
{
   message(sprintf("--- test_vcf.sampleID.to.clinicalTable"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   checkTrue(all(c("EndophenotypeExplorer", "R6") %in% class(etx)))
   tbl.clinical <- etx$vcfSampleID.to.clinicalTable("SM-CTEMF")
   checkEquals(dim(tbl.clinical), c(1, 36))

   rs10200967.00.samples <- c("1036","11241","11249","11254","11275","11312","11337","11349","11370",
                              "11399","11422","11470","14552","1940","1954","1963","68682","6976",
                              "70594","71689","71736","71739","71749","71767","71779","71789","71801",
                              "71824","71860","71884","71886","71900","71901","71903","71904","71908",
                              "71927","71964","71974","71986","76348","894","MAP46251007","MAP50104134",
                              "MAP93787649","SM-CJEG3","SM-CJEG4","SM-CJEGC","SM-CJEGE","SM-CJEGQ",
                              "SM-CJEHR","SM-CJEIU","SM-CJEJI","SM-CJEKL","SM-CJEKW","SM-CJFKO",
                              "SM-CJFKZ","SM-CJFL8","SM-CJFLL","SM-CJFN9","SM-CJFNJ","SM-CJFOC",
                              "SM-CJGH1","SM-CJGI5","SM-CJGIP","SM-CJGJ2","SM-CJGLF","SM-CJGMT",
                              "SM-CJGNF","SM-CJGNJ","SM-CJIXB","SM-CJIYK","SM-CJIZQ","SM-CJIZS",
                              "SM-CJJ17","SM-CJJ1H","SM-CJJ1K","SM-CJJ38","SM-CJK3F","SM-CJK41",
                              "SM-CJK49","SM-CTDQO","SM-CTDQR","SM-CTDQV","SM-CTDQX","SM-CTDR5",
                              "SM-CTDSV","SM-CTDTL","SM-CTDVA","SM-CTDVG","SM-CTECN","SM-CTEE7",
                              "SM-CTEED","SM-CTEFD","SM-CTEGU","SM-CTEID","SM-CTELW","SM-CTEMS")

   tbl.clinical <- etx$vcfSampleID.to.clinicalTable(rs10200967.00.samples)
      # todo: 98 samples ids, only 53 have clinical values
   checkEquals(dim(tbl.clinical), c(53, 36))

} # test_vcf.sampleID.to.clinicalTable
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()

