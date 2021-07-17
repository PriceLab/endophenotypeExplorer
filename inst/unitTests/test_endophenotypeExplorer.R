library(RUnit)
library(EndophenotypeExplorer)
library(plyr)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_readRemoteVCF()
    test_locsToRSID()
    test_getGenoMatrix()
    test_mapSampleIdToPatientAndCohort()
    test_getPatientTables()
    test_vcf.sampleID.to.clinicalTable()
    test_getAggregatedAlleleFrequencies()
    test_gwasLociFrequencies()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   checkTrue(all(c("EndophenotypeExplorer", "R6") %in% class(etx)))
   expected <- "https://igv-data.systemsbiology.net/static/ampad/NIA-1898/chr2.vcf.gz"
   checkEquals(etx$getVcfUrl(), expected)

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_readRemoteVCF <- function(verbose=FALSE)
{
   message(sprintf("--- test_readRemoteVCF"))

   require(VariantAnnotation)
   if(verbose) message(sprintf("packageVersion('VariantAnnotation': %s')", packageVersion("VariantAnnotation")))

   roi <- GRanges(seqnames="2", IRanges(start=127084188, end=127084203))
   url <- "https://igv-data.systemsbiology.net/static/ampad/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_recalibrated_variants/chr2.vcf.gz"
   url <- "https://igv-data.systemsbiology.net/static/ampad/NIA-1898/chr2.vcf.gz"
   x <- readVcf(url, "hg19", roi)

   if(verbose) message(sprintf("size: %d", length(x)))
   checkEquals(dim(geno(x)$GT), c(2, 1894))

} # test_readRemoteVCF
#----------------------------------------------------------------------------------------------------
test_getGenoMatrix <- function()
{
   message(sprintf("--- test_getGenoMatrix"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
      # chr2:127,084,188-127,084,203: 16 base pairs around rs10200967 in bin1 5' utr
   mtx.geno <- etx$getGenoMatrix("2", 127084188, 127094203)
   checkEquals(dim(mtx.geno), c(220, 1894))
   checkTrue(all(c("2:127084193_G/A", "2:127084194_A/C") %in% rownames(mtx.geno)))


} # test_getGenoMatrix
#----------------------------------------------------------------------------------------------------
test_locsToRSID <- function()
{
   message(sprintf("--- test_locsToRSID"))
   list.locs <- get(load(system.file(package="EndophenotypeExplorer", "extdata",
                                     "chromLocs.hg19.for.rsidMapping.RData")))
   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   rsids <- etx$locsToRSID(list.locs, "hg19")
   checkEquals(length(list.locs), length(rsids))
   checkEquals(length(grep(":", rsids)), 69)

     # though the locs are actually hg19, we should be able to make a useless hg38 conversion
   rsids <- etx$locsToRSID(list.locs, "hg38")
   checkEquals(length(list.locs), length(rsids))
   checkEquals(length(grep(":", rsids)), 181)


} # test_locsToRSID
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
test_getGWASTables <- function()
{
   message(sprintf("--- test_vcf.sampleID.to.clinicalTable"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   x <- etx$getGWASTables()
   checkEquals(sort(names(x)), c("gwas.38", "gwas.38.assoc"))
   checkEquals(dim(x$gwas.38), c(38, 13))
   checkEquals(dim(x$gwas.38.assoc), c(3828, 12))
   tbl.bin1.assoc <- subset(x$gwas.38.assoc, gene=="BIN1")
   checkEquals(dim(tbl.bin1.assoc), c(30, 12))

} # test_getGWASTables
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
test_getAggregatedAlleleFrequencies <- function()
{

    message(sprintf("--- test_getAggregatedAlleleFrequencies"))

        # these are all associated with the BIN1 GWAS landmark snp according to
        # Posthuma Genome-wide meta-analysis identifies new loci and functional pathways
        # influencing  Alzheimer's disease risk
        #   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6836675/

    rsids <- c("rs114469753",
               "rs11674295",
               "rs11674920",
               "rs11675014",
               "rs11680116",
               "rs11680117",
               "rs12472671",
               "rs12475915",
               "rs12477427",
               "rs12617835",
               "rs12618593",
               "rs13019485",
               "rs202224448",
               "rs2118506",
               "rs2118508",
               "rs34558561",
               "rs34997561",
               "rs35103166",
               "rs35114168",
               "rs35155471",
               "rs35400498",
               "rs35486865",
               "rs35783077",
               "rs372089992",
               "rs3943703",
               "rs4663105",
               "rs60447541",
               "rs66837244",
               "rs6710467",
               "rs7561528")

    set.seed(17)
       # rs184384746 is rare in Europeans, worhth a look
    rsids <- c("rs184384746", sort(sample(rsids, size=5)))
    tbls <- list()

    etx <- EndophenotypeExplorer$new("BIN1", "hg19")

    for(rsid in rsids){
       printf("--- %s", rsid)
       tbl <- etx$getAggregatedAlleleFrequencies(rsid)
       tbls[[rsid]] <- tbl
       Sys.sleep(1)
       }

   tbl <- do.call(rbind.fill, tbls)
   checkEquals(dim(tbl), c(72, 13))

      # an idiosyncratic test
   tbl.af.rare <- subset(tbl, population=="African American" & (A.freq < 1 | G.freq <1 | C.freq <1 | T.freq < 1))
   checkEquals(tbl.af.rare$rsid, c("rs184384746", "rs114469753"))
      #           rsid ref       population    T    C total     T.freq   C.freq  min.freq  A  G A.freq G.freq
      # 4  rs184384746   C African American    4 2828  2832  0.1412429 99.85876 0.1412429 NA NA     NA     NA
      # 16 rs114469753   T African American 2851    5  2856 99.8249300  0.17507 0.1750700 NA NA     NA     NA

} # test_getAggregatedAlleleFrequencies
#----------------------------------------------------------------------------------------------------
test_gwasLociFrequencies <- function()
{
   message(sprintf("--- test_gwasLociFrequencies"))

   tbl.williams <- get(load(system.file(package="EndophenotypeExplorer", "extdata", "gwas",
                                         "williams-natureNeuroscience2020.RData")))
   dim(tbl.williams)
   checkEquals(dim(tbl.williams), c(58, 2))
   rsids <- unique(grep("^rs", unlist(strsplit(tbl.williams$rsid, ",")), value=TRUE))

   checkEquals(length(unique(rsids)), 57)
   set.seed(17)
   rsids <- rsids[sample(seq_len(57), 5)]
   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   tbls <- list()
   for(rsid in rsids){
      message(sprintf("--- %s", rsid))
      tbl <- etx$getAggregatedAlleleFrequencies(rsid)
      tbls[[rsid]] <- tbl
      Sys.sleep(1)
      }
   tbl <- do.call(rbind.fill, tbls)
   dim(tbl)
   fivenum(tbl$min.freq)
   tbl.rareEuropean <- subset(tbl, population=="European" & min.freq < 1)
   tbl.genesRareEuropean <- subset(tbl.williams, rsid %in% tbl.rareEuropean$rsid)
   checkEquals("CNTNAP2", tbl.genesRareEuropean$locusOrGene)

} # test_gwasLociFrequencies
#----------------------------------------------------------------------------------------------------

if(!interactive())
   runTests()
