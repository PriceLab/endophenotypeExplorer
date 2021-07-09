library(RUnit)
library(EndophenotypeExplorer)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_getGenoMatrix()
    test_mapSampleIdToPatientAndCohort()
    test_getPatientTables()
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
test_getGenoMatrix <- function()
{
   message(sprintf("--- test_getGenoMatrix"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
      # chr2:127,084,188-127,084,203: 16 base pairs around rs10200967 in bin1 5' utr
   mtx.geno <- etx$getGenoMatrix("2", 127084188, 127084203)
   checkEquals(dim(mtx.geno), c(2, 1894))
   checkEquals(rownames(mtx.geno), c("2:127084193_G/A", "2:127084194_A/C"))

} # test_getGenoMatrix
#----------------------------------------------------------------------------------------------------
test_mapSampleIdToPatientAndCohort <- function()
{
   message(sprintf("--- test_mapSampleIdToPatientAndCohort"))
   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   tbl.map <- etx$getIdMap()
   checkEquals(dim(tbl.map), c(1845, 3))
   checkEquals(colnames(tbl.map), c("vcf", "patient", "cohort"))

   map <- as.list(etx$mapSampleIdToPatientAndCohort("SM-CJGH1"))
   checkEquals(map, list(vcf="SM-CJGH1", patient="R7025378", cohort="rosmap"))

} # test_mapSampleIdToPatientAndCohort
#----------------------------------------------------------------------------------------------------
test_vcf.sampleID.to.clinicalTable <- function()
{
   message(sprintf("--- test_vcf.sampleID.to.clinicalTable"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   checkTrue(all(c("EndophenotypeExplorer", "R6") %in% class(etx)))

      # choose just one
   tbl.clinical <- etx$vcfSampleID.to.clinicalTable("SM-CTEMF")
   checkEquals(dim(tbl.clinical), c(1, 18))
   checkEquals(colnames(tbl.clinical), c("projid", "Study", "msex", "educ", "race", "spanish",
                                         "apoe_genotype", "age_at_visit_max", "age_first_ad_dx",
                                         "age_death", "cts_mmse30_first_ad_dx", "cts_mmse30_lv",
                                         "pmi", "braaksc", "ceradsc", "cogdx", "dcfdx_lv", "individualID"))
      # spot check a few important fields
   checkEquals(tbl.clinical$apoe_genotype, 33)
   checkEquals(tbl.clinical$braaksc, 4)
   checkEquals(tbl.clinical$ceradsc, 2)

   problem.sampleID <- "71729"
   tbl.clinical <- etx$vcfSampleID.to.clinicalTable(problem.sampleID)
   checkEquals(nrow(tbl.clinical), 0)

} # test_vcf.sampleID.to.clinicalTable
#----------------------------------------------------------------------------------------------------
test_getPatientTables <- function()
{
   message(sprintf("--- test_getPatientTables"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")

      #--------------------------
      # first, get full tables
      #--------------------------

   tbl.pt.sinai <- etx$get.sinai.patient.table(NA)
   checkEquals(dim(tbl.pt.sinai), c(377, 20))

   tbl.pt.rosmap <- etx$get.rosmap.patient.table(NA)
   checkEquals(dim(tbl.pt.rosmap), c(3584, 18))

   tbl.pt.mayo <- etx$get.mayo.patient.table(NA)
   checkEquals(dim(tbl.pt.mayo), c(370, 19))

      #---------------------------------
      # now just one patient at a time
      #---------------------------------

   tbl.pt.sinai <- etx$get.sinai.patient.table("AMPAD_MSSM_0000026992")
   checkEquals(dim(tbl.pt.sinai), c(1, 20))

   tbl.pt.rosmap <- etx$get.rosmap.patient.table("R1977848")
   checkEquals(dim(tbl.pt.rosmap), c(1, 18))

   tbl.pt.mayo <- etx$get.mayo.patient.table(1000)
   checkEquals(dim(tbl.pt.mayo), c(1, 19))

} # test_getPatientTables
#----------------------------------------------------------------------------------------------------
test_vcf.sampleID.to.clinicalTable <- function()
{
   message(sprintf("--- test_vcf.sampleID.to.clinicalTable"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   tbl.map <- etx$getIdMap()

   vcf.ids <- c(subset(tbl.map, cohort=="sinai")$vcf[1],
                subset(tbl.map, cohort=="rosmap")$vcf[1],
                subset(tbl.map, cohort=="mayo")$vcf[1])

   tbls <- lapply(vcf.ids, function(id) etx$vcfSampleID.to.clinicalTable(id))
   tbl <- do.call(rbind, tbls)
   checkEquals(dim(tbl), c(3, 12))
   checkTrue(is.numeric(tbl$ageAtDeath))

   # as.data.frame(t(tbl))
   set.seed(17)
   vcf.ids <- sample(tbl.map$vcf, size=100)
   tbls <- lapply(vcf.ids, function(id) etx$vcfSampleID.to.clinicalTable(id))
   tbl.all <- do.call(rbind, tbls)
   checkEquals(dim(tbl.all), c(100, 12))

   checkEqualsNumeric(mean(tbl.all$braak, na.rm=TRUE), 3.9, tolerance=0.2)
   checkEqualsNumeric(mean(tbl.all$pmi, na.rm=TRUE),    78, tolerance=0.2)
   checkEqualsNumeric(mean(tbl.all$cogdx, na.rm=TRUE), 2.5, tolerance=0.2)
   checkEqualsNumeric(mean(tbl.all$cerad, na.rm=TRUE), 2.3, tolerance=0.2)

   checkEqualsNumeric(fivenum(tbl.all$ageAtDeath),
                      c(66.00000, 84.92300, 89.97878, 90.00000, 90.00000),
                      tolerance=0.5)

   checkEquals(length(which(is.na(tbl.all$cogdx))), 32)

} # test_vcf.sampleID.to.clinicalTable
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
