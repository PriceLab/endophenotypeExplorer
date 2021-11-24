library(RUnit)
library(EndophenotypeExplorer)
library(plyr)
# source("~/github/endophenotypeExplorer/R/getExpressionMatrices.R")
#----------------------------------------------------------------------------------------------------
init.snpLocs <- function() x <- EndophenotypeExplorer$new(NA, NA, initialize.snpLocs=TRUE)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_getExpressionMatrix()
    test_standardizeRosmapPatientTable ()
    test_standardizeSinaiPatientTable()
    test_standardizeMayoPatientTable()

    test_sinaiMapping()
    test_rosmapMapping()

    test_getClinicalTable()
    test_readRemoteVCF()
    test_locsToRSID()
    test_rsidToLocs()
    test_getGenoMatrix()
    test_mapSampleIdToPatientAndCohort()
    test_getPatientTables()
    test_sampleID.to.clinicalTable()
    test_getAggregatedAlleleFrequencies()
    test_gwasLociFrequencies()
    test_gtexTissueExpression()
    test_getEQTLsForGene()

    test_splitExpressionMatrixByMutationStatusAtRSID_mayo()
    test_splitExpressionMatrixByMutationStatusAtRSID_sinai()
    test_splitExpressionMatrixByMutationStatusAtRSID_rosmap()

    test_subsetAndRelabelGenoMatrixByPatientIDs()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
   message(sprintf("--- test_ctor"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   checkTrue(all(c("EndophenotypeExplorer", "R6") %in% class(etx)))
   expected <- "https://igv-data.systemsbiology.net/static/ampad/NIA-1898/chr2.vcf.gz"
   checkEquals(etx$getVcfUrl(), expected)

   etx <- EndophenotypeExplorer$new(NA, NA)
   etx$setTargetGene("BIN1", "hg19")
   checkTrue(all(c("EndophenotypeExplorer", "R6") %in% class(etx)))
   checkEquals(etx$getVcfUrl(), expected)

   etx$setTargetGene("BACH1", "hg38")
   expected <- "https://igv-data.systemsbiology.net/static/ampad/NIA-1898/chr21.vcf.gz"
   checkEquals(etx$getVcfUrl(), expected)


} # test_ctor
#----------------------------------------------------------------------------------------------------
test_getExpressionMatrix <- function()
{
    message(sprintf("--- test_getExpressionMatrix"))

    etx <- EndophenotypeExplorer$new("BIN1", "hg38")
    codes <- names(etx$get.rna.matrix.codes())
    checkTrue(length(codes) >= 10)
      # check a few
    checkTrue(all(c("max-tcx", "sage-eqtl-cer", "old-rosmap", "max-rosmap") %in%  codes))

      # now loop through all, get the matrices, check for reasonable dimensions
    for(code in codes){
        mtx <- etx$get.rna.matrix(code)
        printf("--- %s: %d x %d", code, nrow(mtx), ncol(mtx))
        checkTrue(nrow(mtx) > 14000)
        checkTrue(ncol(mtx) > 250)
        }

} # test_getExpressionMatrix
#----------------------------------------------------------------------------------------------------
test_standardizeMayoPatientTable <- function()
{
    message(sprintf("--- test_standardizeMayoPatientTable"))

    etx <- EndophenotypeExplorer$new("BIN1", "hg19")

    dir <- system.file(package="EndophenotypeExplorer", "extdata", "clinical")
    tbl.clinical.mayo <- read.table(file.path(dir, "MayoRNAseq_individual_metadata.csv"),
                                    sep=",", header=TRUE, as.is=TRUE)
    tbl.std <- etx$standardizeMayoPatientTable(tbl.clinical.mayo)
    checkEquals(nrow(tbl.std), nrow(tbl.clinical.mayo))

    checkEquals(colnames(tbl.std), etx$getStandardClinicalColumnNames())

} # test_standardizeMayoPatientTable
#----------------------------------------------------------------------------------------------------
test_standardizeRosmapPatientTable <- function()
{
    message(sprintf("--- test_standardizeRosmapPatientTable"))

    etx <- EndophenotypeExplorer$new("BIN1", "hg19")

    dir <- system.file(package="EndophenotypeExplorer", "extdata", "clinical")
    file.exists(dir)
    full.path <- file.path(dir, "ROSMAP_clinical.csv")
    file.exists(full.path)
    tbl.clinical.rosmap <- read.table(full.path, sep=",", header=TRUE, as.is=TRUE)
    checkEquals(dim(tbl.clinical.rosmap), c(3583, 18))

    tbl.std <- etx$standardizeRosmapPatientTable(tbl.clinical.rosmap)
    checkEquals(nrow(tbl.std), nrow(tbl.clinical.rosmap))
    checkEquals(colnames(tbl.std), etx$getStandardClinicalColumnNames())

} # test_standardizeRosmapPatientTable
#----------------------------------------------------------------------------------------------------
test_standardizeSinaiPatientTable <- function()
{
    message(sprintf("--- test_standardizeSinaiPatientTable"))

    etx <- EndophenotypeExplorer$new("BIN1", "hg19")

    dir <- system.file(package="EndophenotypeExplorer", "extdata", "clinical")
    tbl.clinical.sinai <- read.table(file.path(dir, "MSBB_individual_metadata.csv"),
                                    sep=",", header=TRUE, as.is=TRUE)
    checkEquals(dim(tbl.clinical.sinai), c(377, 20))
    tbl.std <- etx$standardizeSinaiPatientTable(tbl.clinical.sinai)
    checkEquals(nrow(tbl.std), nrow(tbl.clinical.sinai))

    checkEquals(colnames(tbl.std), etx$getStandardClinicalColumnNames())

} # test_standardizeSinaiPatientTable
#----------------------------------------------------------------------------------------------------
# a successful sinai mapping has a vcf and an rna-seq sample id for every sinai patient
test_sinaiMapping <- function()
{
    message(sprintf("--- test_sinaiMapping"))

    dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
    file.sinai <- "mtx.mssm.rnaseq-residual-eqtl-geneSymbols-patients-16346x753.RData"
    mtx.sinai <- get(load(file.path(dir, file.sinai)))
    checkEquals(dim(mtx.sinai), c(16346, 753))

    etx <- EndophenotypeExplorer$new("PTK2B", "hg38")
    tbl.map <- etx$getIdMap()
    pts.sinai <- unique(subset(tbl.map, study=="sinai" & !is.na(patient))$patient)
    checkEquals(length(pts.sinai), 358)
    checkTrue(all(grepl("AMPAD_MSSM_", pts.sinai)))

    rna.samples <- unique(subset(tbl.map, study=="sinai" & assay=="rnaseq")$sample)
    checkEquals(length(rna.samples), 1026)
    rna.ids.hB <- grep("^hB_RNA", rna.samples, value=TRUE)
    rna.ids.BM <- grep("^BM_",    rna.samples, value=TRUE)
    checkEquals(length(rna.ids.hB), 597)
    checkEquals(length(rna.ids.BM), 429)
    checkEquals(length(rna.ids.hB) + length(rna.ids.BM), length(rna.samples))

    vcf.samples <- unique(subset(tbl.map, study=="sinai" & assay=="vcf")$sample)
    length(vcf.samples)
    checkEquals(length(vcf.samples), 345)

       #--------------------------------------------------
       # an important query:
       #   obtain the rnaseq sample ids corresponding to
       #   vcf sample ids, mapped through the patient
       #--------------------------------------------------

    vcf.patients <- unique(subset(tbl.map, study=="sinai" & assay=="vcf")$patient)
    checkEquals(length(vcf.patients), 341)

    rna.patients <- unique(subset(tbl.map, study=="sinai" & assay=="rnaseq")$patient)
    checkEquals(length(rna.patients), 302)

    both.patients <- intersect(vcf.patients, rna.patients)
    length(both.patients)  # 284

       # a typical use of this table:
       #  given a vcf sample ids, determine the patients, return
       #  the corresponding rnaseq, identify the (likely) subset of
       #  patients with both kinds of data

    set.seed(17)
    poi <- sort(rna.patients[sample(seq_len(length(rna.patients)), size=30)])
    vcf.samples <- subset(tbl.map, patient %in% poi & study=="sinai" & assay=="vcf")$sample

      # find the patients
    patients <- unique(subset(tbl.map, sample %in% vcf.samples & study=="sinai")$patient)
    checkEquals(length(patients), 29)

      # now pick out the columns of the rnaseq matrix, which have already
      # been mapped to patient
    mtx.rnaseq.colnames <- intersect(patients, colnames(mtx.sinai))
    checkEquals(length(mtx.rnaseq.colnames), 23)

} # test_sinaiMapping
#----------------------------------------------------------------------------------------------------
test_rosmapMapping <- function()
{
    message(sprintf("--- test_rosmapMapping"))

    dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
    file.rosmap <- "mtx.rosmap.rnaseq-residual-eqtl-geneSymbols-patients-15582x632.RData"
    mtx.rosmap <- get(load(file.path(dir, file.rosmap)))
    checkEquals(dim(mtx.rosmap), c(15582, 632))

    etx <- EndophenotypeExplorer$new("PTK2B", "hg38")
    tbl.map <- etx$getIdMap()
    pts.rosmap <- unique(subset(tbl.map, study=="rosmap" & !is.na(patient))$patient)
    checkEquals(length(pts.rosmap), 1223)
    checkTrue(all(grepl("^R", pts.rosmap)))

    rna.samples <- unique(subset(tbl.map, study=="rosmap" & assay=="rnaseq")$sample)
    checkEquals(length(rna.samples), 632)
       # rna.sample ids look like this
       # "01_120405" "02_120405" "03_120405" "04_120405" "05_120405"
       # where the first tokens appear to range from "01" to "958"
       # plus "R24" and "redo4"
    # the second tokens cluster like this:
    # head(table(unlist(lapply(strsplit((rna.samples), "_"), "[", 2))))
    #
    # 120405 120410 120411 120416 120417 120418
    #  5      6      8     17     34     16
    # i'm not sure what these group number refer to

    most.populous.group <- names(sort(table(unlist(lapply(
        strsplit((rna.samples), "_"), "[", 2))), decreasing=TRUE))[1]
    checkEquals(most.populous.group, "120430") # 37 members

    vcf.samples <- unique(subset(tbl.map, study=="rosmap" & assay=="vcf")$sample)
    checkEquals(length(vcf.samples), 1151)

       #--------------------------------------------------
       # an important query:
       #   obtain the rnaseq sample ids corresponding to
       #   vcf sample ids, mapped through the patient
       #--------------------------------------------------

    vcf.patients <- unique(subset(tbl.map, study=="rosmap" & assay=="vcf")$patient)
    checkEquals(length(vcf.patients), 1144)

    rna.patients <- unique(subset(tbl.map, study=="rosmap" & assay=="rnaseq")$patient)
    checkEquals(length(rna.patients), 632)

    both.patients <- intersect(vcf.patients, rna.patients)
    checkEquals(length(both.patients), 553)

       # a typical use of this table:
       #   given vcf sample ids, determine the patients, return
       #   the corresponding rnaseq, identify the (likely) subset of
       #   patients with both kinds of data

    set.seed(17)
    poi <- sort(rna.patients[sample(seq_len(length(rna.patients)), size=30)])
    checkEquals(length(poi), 30)

       # here are our vcf sample ids to try out:
    vcf.samples <- subset(tbl.map, patient %in% poi &
                                   study=="rosmap" &
                                   assay=="vcf")$sample
    checkEquals(length(vcf.samples), 27)

      # find the patients
    patients <- unique(subset(tbl.map, sample %in% vcf.samples & study=="rosmap")$patient)
    checkEquals(length(patients), 27)

      # now pick out the columns of the rnaseq matrix, which have already
      # been mapped to patient
    mtx.rnaseq.colnames <- intersect(patients, colnames(mtx.rosmap))
    checkEquals(length(mtx.rnaseq.colnames), 27)

} # test_rosmapMapping
#----------------------------------------------------------------------------------------------------
test_getClinicalTable <- function()
{
   message(sprintf("--- test_getClinicalTable"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19",
                                    defer.setupClinicalData.toSupportTesting=TRUE)
   etx$setupClinicalData()

   tbl <- etx$getClinicalTable()
   checkEquals(colnames(tbl), etx$getStandardClinicalColumnNames())

   studies <- as.list(table(tbl$study))
   checkEquals(studies$ROS, 1451)
   checkEquals(studies$MAP, 2132)
   checkEquals(studies$MSSM, 377)
   checkEquals(studies$MayoBrainBank, 370)

   cogdx <- as.list(table(tbl$cogdx))
   names(cogdx) <- paste0("cogdx.", names(cogdx))

   checkEquals(cogdx$cogdx.1, 586)
   checkEquals(cogdx$cogdx.2, 404)
   checkEquals(cogdx$cogdx.3, 33)
   checkEquals(cogdx$cogdx.4, 674)
   checkEquals(cogdx$cogdx.5, 94)
   checkEquals(cogdx$cogdx.6, 30)

} # test_getClinicalTable
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
   loc.start <- 127084188
   loc.end   <- 127094203
   x <- system.time(mtx.geno <- etx$getGenoMatrix("2", loc.start, loc.end))
   checkEquals(dim(mtx.geno), c(220, 1894))
   printf("--- retrieved %d variants, %d samples, in %d seconds across %d kb bases on chr2",
          nrow(mtx.geno), ncol(mtx.geno),
          as.integer(x[["elapsed"]]),
          as.integer((loc.end - loc.start)/1000))
   checkTrue(all(c("2:127084193_G/A", "2:127084194_A/C") %in% rownames(mtx.geno)))


     #--------------------------------------------------------------------------------
     # divide the region into two halves, should get the same matrix back
     #--------------------------------------------------------------------------------

   mtx.geno.2 <- etx$getGenoMatrix(chrom=c("2", "2"),
                                   start=c(127084188, 127084196),
                                   end  =c(127084195, 127094203))
   checkEquals(dim(mtx.geno.2), dim(mtx.geno))

     #--------------------------------------------------------------------------------
     # query across a much larger area, that which spans all the rosmap/ampad eQTLs
     # for NDUFS2
     #--------------------------------------------------------------------------------
   etx$setTargetGene("NDUFS2", "hg19")
   tbl.eqtls.ndufs2 <- etx$getEQTLsForGene()
   checkEquals(dim(tbl.eqtls.ndufs2), c(11310, 10))
   tbl.eqtls.sub <- subset(tbl.eqtls.ndufs2, study=="ampad-rosmap" & pvalue < 1e-5)
   dim(tbl.eqtls.sub)


   loc.start <- min(tbl.eqtls.sub$hg19)
   loc.end  <-  max(tbl.eqtls.sub$hg19)
   loc.chrom  <- "1"

   x <- system.time(mtx.geno <- etx$getGenoMatrix(loc.chrom, loc.start, loc.end))
   checkEquals(dim(mtx.geno), c(2155, 1894))
   printf("--- retrieved %d variants, %d samples, in %d seconds across %d kb bases on chr2",
          nrow(mtx.geno), ncol(mtx.geno),
          as.integer(x[["elapsed"]]),
          as.integer((loc.end - loc.start)/1000))




     #--------------------------------------------------------------------------------
     # now, on a different chromosome, query by rsid, then write rsids back into the mtx rownames
     #--------------------------------------------------------------------------------

   etx <- EndophenotypeExplorer$new("NDUFS2", "hg19")
   rsids <- c("rs11576415", "rs11584174", "rs12753774", "rs12754503")
   x <- system.time(mtx.geno <- etx$getGenoMatrixByRSID(rsids))
   printf("--- retrieved %d variants, %d samples, in %d seconds",
          nrow(mtx.geno), ncol(mtx.geno),
          as.integer(x[["elapsed"]]))



   new.names <- etx$locsToRSID(rownames(mtx.geno), "hg19")
   checkEquals(names(new.names), rownames(mtx.geno))
        # the locations were sorted before the matrix was retreived.
        # but the rsids are not in location order
        # thus the first check will fail, and so is negated in order to pass
   checkTrue(!all(as.character(new.names) == rsids))
   checkEquals(sort(as.character(new.names)), sort(rsids))
   rownames(mtx.geno) <- as.character(new.names)
   checkEquals(sort(rownames(mtx.geno)), sort(rsids))

} # test_getGenoMatrix
#----------------------------------------------------------------------------------------------------
test_getGenoMatrixByRSID_atScale <- function()
{
    message(sprintf("--- test_getGenoMatrixByRSID_atScale"))
    etx$setTargetGene("NDUFS2", "hg19")
    tbl.eqtls <- etx$getEQTLsForGene()

    time.001 <- system.time(mtx.geno.001 <- etx$getGenoMatrixByRSID(tbl.eqtls$rsid[1]))
    time.010 <- system.time(mtx.geno.001 <- etx$getGenoMatrixByRSID(tbl.eqtls$rsid[1:10]))
    time.100 <- system.time(mtx.geno.001 <- etx$getGenoMatrixByRSID(tbl.eqtls$rsid[1:100]))

    rsids.100 <- tbl.eqtls$rsid[1:100]
    library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
    library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
    time.100.snpsById.hg19 <- system.time(x <- lapply(rsids.100, function(rsid)
                snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, rsid, ifnotfound="drop")))
    time.100.snpsById.hg38 <- system.time(x <- lapply(rsids.100, function(rsid)
                snpsById(SNPlocs.Hsapiens.dbSNP151.GRCh38, rsid, ifnotfound="drop")))
    time.100.snpsById.hg38.oneCall <- system.time(
                snpsById(SNPlocs.Hsapiens.dbSNP151.GRCh38, rsids.100, ifnotfound="drop"))

    system.time(etx$rsidToLoc(rsids.100))

} # test_getGenoMatrixByRSID_atScale
#----------------------------------------------------------------------------------------------------
test_getGenoMatrix_multipleGenes <- function()
{
   message(sprintf("--- test_getGenoMatrix_multipleGenes"))

   etx <- EndophenotypeExplorer$new(NA, NA, initialize.snpLocs=TRUE)
   tbl.gwas <- etx$getGWASTables()$gwas.38

   genes <- sort(tbl.gwas$geneSymbol[c(1,38)])
   checkEquals(genes, c("AGRN", "APP"))
      # make sure these are different chromosomes to provide a good test
   checkTrue(tbl.gwas$chrom[1] != tbl.gwas$chrom[38])

   etx$setTargetGene("AGRN", "hg19")
   checkEquals(etx$getVcfUrl(),
               "https://igv-data.systemsbiology.net/static/ampad/NIA-1898/chr1.vcf.gz")
   mtx.geno.1 <- etx$getGenoMatrixByRSID(tbl.gwas$leadVariant[1])
   checkEquals(dim(mtx.geno.1), c(1, 1894))

   etx$setTargetGene("APP", "hg19")
   checkEquals(etx$getVcfUrl(),
               "https://igv-data.systemsbiology.net/static/ampad/NIA-1898/chr21.vcf.gz")
   mtx.geno.38 <- etx$getGenoMatrixByRSID(tbl.gwas$leadVariant[38])
   checkEquals(dim(mtx.geno.38), c(1, 1894))

      # make sure we did not get the same matrix back twice
   checkTrue(as.list(table(mtx.geno.1))[["0/1"]] != as.list(table(mtx.geno.38))[["0/1"]])

} # test_getGenoMatrix_multipleGenes
#----------------------------------------------------------------------------------------------------
test_locsToRSID <- function()
{
   message(sprintf("--- test_locsToRSID"))
   etx <- EndophenotypeExplorer$new("BIN1", "hg19", initialize.snpLocs=TRUE)

     # get a set of rsids with hg19 and hg38 locs
   tbl.eqtl <- etx$getEQTLsForGene()
   tbl.bin1 <-  head(subset(tbl.eqtl, study=="ampad-rosmap"))
   locs.hg19 <- sprintf("%s:%d_X/Y", tbl.bin1$chrom, tbl.bin1$hg19)
   locs.hg19 <- sub("chr", "", locs.hg19)
   locs.hg38 <- sprintf("%s:%d_X/Y", tbl.bin1$chrom, tbl.bin1$hg38)
   locs.hg38 <- sub("chr", "", locs.hg38)
   rsids.orig<- tbl.bin1$rsid

      #----------------------------------------
      # first, hg19 locs
      #----------------------------------------

   rsids <- etx$locsToRSID(locs.hg19, "hg19")
   checkEquals(names(rsids), locs.hg19)
   checkTrue(all(rsids.orig %in% rsids))

      #----------------------------------------
      # now, hg38 locs
      #----------------------------------------

   rsids <- etx$locsToRSID(locs.hg38, "hg38")
   checkEquals(length(locs.hg38), length(rsids))
   checkTrue(all(rsids.orig %in% rsids))


} # test_locsToRSID
#----------------------------------------------------------------------------------------------------
test_rsidToLocs <- function()
{
   message(sprintf("--- test_rsidToLocs"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")

   rsid.list <- "rs4575098"
   tbl.locs <- etx$rsidToLoc(rsid.list)
   checkEquals(tbl.locs$chrom, "1")
   checkEquals(tbl.locs$hg19, 161155392)
   checkEquals(tbl.locs$hg38, 161185602)
   checkEquals(tbl.locs$rsid, "rs4575098")

   rsid.list <- c("rs114360492", "rs4351014", "rs74615166", "rs6504163", "rs9381040")
   system.time(tbl.locs <- etx$rsidToLoc(rsid.list))

   new.order <- match(rsid.list, tbl.locs$rsid)
   tbl.locs <- tbl.locs[new.order,]
   checkEquals(tbl.locs$rsid, rsid.list)
   checkEquals(tbl.locs$chrom, c("7", "4", "15", "17", "6"))
   checkEquals(tbl.locs$hg19,  c(145950029, 11027619, 64725490, 61545779, 41154650))
   checkEquals(tbl.locs$hg38,  c(146252937, 11025995, 64433291, 63468418, 41186912))

   rsid.list <- c("rs4575098", "bogus")
   tbl.locs <- etx$rsidToLoc(rsid.list)
   checkEquals(dim(tbl.locs), c(1, 4))

   rsid.list <- c("rs4575099", "bogus")   # rs4575099 does not exist
   tbl.locs <- etx$rsidToLoc(rsid.list)
   checkEquals(dim(tbl.locs), c(0, 4))

} #  test_rsidToLocs()
#----------------------------------------------------------------------------------------------------
test_mapSampleIdToPatientAndCohort <- function()
{
   message(sprintf("--- test_mapSampleIdToPatientAndCohort"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   tbl.map <- etx$getIdMap()
   checkEquals(dim(tbl.map), c(4026, 4))
   checkEquals(colnames(tbl.map), c("sample", "patient", "study", "assay"))

   x <- etx$mapSampleIdToPatientAndCohort("SM-CJGH1")
   checkEquals(dim(x), c(1, 4))
   checkEquals(x$sample,  "SM-CJGH1")
   checkEquals(x$patient, "R7025378")
   checkEquals(x$study,   "rosmap")
   checkEquals(x$assay,    "vcf")

} # test_mapSampleIdToPatientAndCohort
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
   checkEquals(dim(tbl.pt.rosmap), c(3583, 18))

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
   message(sprintf("--- test_getGWASTables"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   x <- etx$getGWASTables()
   checkEquals(sort(names(x)), c("gwas.38", "gwas.38.assoc"))
   checkEquals(dim(x$gwas.38), c(38, 13))
   checkEquals(dim(x$gwas.38.assoc), c(3828, 12))
   tbl.bin1.assoc <- subset(x$gwas.38.assoc, gene=="BIN1")
   checkEquals(dim(tbl.bin1.assoc), c(30, 12))

} # test_getGWASTables
#----------------------------------------------------------------------------------------------------
test_sampleID.to.clinicalTable <- function()
{
   message(sprintf("--- test_sampleID.to.clinicalTable"))

   etx <- EndophenotypeExplorer$new("BIN1", "hg19")
   tbl.map <- etx$getIdMap()
   tbl.clinical <- etx$getClinicalTable()

   sample.ids <- c(subset(tbl.map, study=="sinai")$sample[1],
                   subset(tbl.map, study=="rosmap")$sample[1],
                   subset(tbl.map, study=="mayo")$sample[1])

   tbls <- lapply(sample.ids, function(id) etx$sampleID.to.clinicalTable(id))
   tbl <- do.call(rbind, tbls)
   checkEquals(dim(tbl), c(3, 12))
   checkTrue(is.numeric(tbl$ageAtDeath))

   set.seed(17)
   sample.ids <- sample(tbl.map$sample, size=100)
   tbls <- lapply(sample.ids, function(id) etx$sampleID.to.clinicalTable(id))
   tbl.all <- do.call(rbind, tbls)
   checkEquals(dim(tbl.all), c(100, 12))

   checkEqualsNumeric(mean(tbl.all$braak, na.rm=TRUE), 3.9, tolerance=0.2)
   checkEqualsNumeric(mean(tbl.all$pmi, na.rm=TRUE),   122.9, tolerance=0.2)
   checkEqualsNumeric(mean(tbl.all$cogdx, na.rm=TRUE), 2.5, tolerance=0.2)
   checkEqualsNumeric(mean(tbl.all$cerad, na.rm=TRUE), 2.3, tolerance=0.2)

   checkEqualsNumeric(fivenum(tbl.all$ageAtDeath),
                      c(66.00000, 84.92300, 89.97878, 90.00000, 90.00000),
                      tolerance=0.5)

   checkEquals(length(which(is.na(tbl.all$cogdx))), 52)

} # test_sampleID.to.clinicalTable
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
    length(rsids)
    tbls <- list()

    etx <- EndophenotypeExplorer$new("PTK2B", "hg38")

    for(rsid in rsids){
       printf("--- %s", rsid)
       tbl <- etx$getAggregatedAlleleFrequencies(rsid, quiet=FALSE)
       tbls[[rsid]] <- tbl
       Sys.sleep(2)
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
      Sys.sleep(2)
      }
   tbl <- do.call(rbind.fill, tbls)
   dim(tbl)
   fivenum(tbl$min.freq)
   tbl.rareEuropean <- subset(tbl, population=="European" & min.freq < 1)
   tbl.genesRareEuropean <- subset(tbl.williams, rsid %in% tbl.rareEuropean$rsid)
   checkEquals("CNTNAP2", tbl.genesRareEuropean$locusOrGene)

} # test_gwasLociFrequencies
#----------------------------------------------------------------------------------------------------
test_gtexTissueExpression <- function()
{
    message(sprintf("--- test_gtexTissueExpression"))
    targetGene <- "BIN1"
    etx <- EndophenotypeExplorer$new(targetGene, "hg19")
    tbl.gtex <- etx$getTissueExpression(targetGene, "brain")
    checkEquals(dim(tbl.gtex), c(13, 6))
    checkEquals(colnames(tbl.gtex),
                c("datasetId", "gencodeId", "geneSymbol", "median", "tissueSiteDetailId", "unit"))
    checkTrue(all(tbl.gtex$geneSymbol == targetGene))
    checkEqualsNumeric(mean(tbl.gtex$median), 140, tolerance=1)

} # test_gtexTissueExpression
#----------------------------------------------------------------------------------------------------
test_getEQTLsForGene <- function()
{
    message(sprintf("--- test_getEQTLsForGene"))

    targetGene <- "NDUFS2"
    etx <- EndophenotypeExplorer$new(targetGene, "hg19", initialize.snpLocs=FALSE)
    tbl.eQTL <- etx$getEQTLsForGene()
    checkTrue(nrow(tbl.eQTL) > 10000)
    checkTrue(ncol(tbl.eQTL) >= 10)

    min.pval <- min(tbl.eQTL$pvalue)
    checkTrue("rs1136224" %in% subset(tbl.eQTL, pvalue == min.pval)$rsid)

       # has the table been sorted?
    checkEquals(tbl.eQTL$rsid[1], "rs1136224")

    tbl.sig <- subset(tbl.eQTL, pvalue <= 0.01)
    dim(tbl.sig)
    checkTrue(nrow(tbl.sig) < 400 & nrow(tbl.sig) > 300)

} # test_getEQTLsForGene
#----------------------------------------------------------------------------------------------------
test_splitExpressionMatrixByMutationStatusAtRSID_mayo <- function()
{
   message(sprintf("--- test_splitExpressionMatrixByMutationStatusAtRSID_mayo"))

   dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
   file.cer <- "mtx.mayo.cer.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x255.RData"
   mtx.cer <- get(load(file.path(dir, file.cer)))

   targetGene <- "PTK2B"
   rsid <- "rs28834970"

   etx <- EndophenotypeExplorer$new(targetGene, "hg38")

       #------------------------------------------------------------
       # split mayo cerebellum (mtx.cer) into 4 matrices by genotype
       #------------------------------------------------------------

   x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx.cer, rsid, study.name="mayo")
   checkEquals(sort(names(x)), c("genotypes.rna", "genotypes.vcf",
                                 "het", "hom", "mut", "wt"))
   checkEquals(x$genotypes.rna, list(wt=95, mut=160, het=124, hom=36))
   checkEquals(x$genotypes.vcf, list(wt=129, mut=220, het=167, hom=53))

   checkEquals(sum(ncol(x$wt), ncol(x$mut)), ncol(mtx.cer))
   checkEquals(sum(ncol(x$het), ncol(x$hom)), ncol(x$mut))

   checkEquals(dim(x$wt), c(17009, 95))
   checkEquals(dim(x$mut), c(17009, 160))
   checkEquals(dim(x$het), c(17009, 124))
   checkEquals(dim(x$hom), c(17009, 36))

       #------------------------------------------------------------------
       # split mayo temporal cortex (mtx.tcx) into 4 matrices by genotype
       #------------------------------------------------------------------

   file.tcx <- "mtx.mayo.tcx.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x257.RData"
   mtx.tcx <- get(load(file.path(dir, file.tcx)))
   checkEquals(dim(mtx.tcx), c(17009, 257))

   x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx.tcx, rsid, study.name="mayo")
   checkEquals(sort(names(x)), c("genotypes.rna", "genotypes.vcf",
                                 "het", "hom", "mut", "wt"))
   checkEquals(x$genotypes.rna, list(wt=94, mut=163, het=122, hom=41))
   checkEquals(x$genotypes.vcf, list(wt=129, mut=220, het=167, hom=53))

   checkEquals(sum(ncol(x$wt), ncol(x$mut)), ncol(mtx.tcx))
   checkEquals(sum(ncol(x$het), ncol(x$hom)), ncol(x$mut))

   checkEquals(dim(x$wt), c(17009, 94))
   checkEquals(dim(x$mut), c(17009, 163))
   checkEquals(dim(x$het), c(17009, 122))
   checkEquals(dim(x$hom), c(17009, 41))

} # test_splitExpressionMatrixByMutationStatusAtRSID
#----------------------------------------------------------------------------------------------------
test_splitExpressionMatrixByMutationStatusAtRSID_sinai <- function()
{
   message(sprintf("--- test_splitExpressionMatrixByMutationStatusAtRSID_sinai"))

   dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
   file.sinai <- "mtx.mssm.rnaseq-residual-eqtl-geneSymbols-patients-16346x753.RData"
   mtx <- get(load(file.path(dir, file.sinai)))
   checkEquals(dim(mtx), c(16346, 753))

   targetGene <- "PTK2B"
   rsid <- "rs28834970"

   etx <- EndophenotypeExplorer$new(targetGene, "hg38")
   x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx, rsid, study.name="sinai")

   checkEquals(sort(names(x)), c("genotypes.rna", "genotypes.vcf",
                                 "het", "hom", "mut", "wt"))
   checkEquals(x$genotypes.rna, list(wt=108, mut=128, het=103, hom=25))
   checkEquals(x$genotypes.vcf, list(wt=160, mut=181, het=143, hom=38))

     # the sinai matrix has lots of as yet poorly understood samples
   checkEquals(ncol(mtx) - sum(ncol(x$wt), ncol(x$mut)), 517)
   checkEquals(sum(ncol(x$het), ncol(x$hom)), ncol(x$mut))

   checkEquals(dim(x$wt), c(16346,  108))
   checkEquals(dim(x$mut), c(16346, 128))
   checkEquals(dim(x$het), c(16346, 103))
   checkEquals(dim(x$hom), c(16346, 25))

} # test_splitExpressionMatrixByMutationStatusAtRSID_sinai
#----------------------------------------------------------------------------------------------------
test_splitExpressionMatrixByMutationStatusAtRSID_rosmap <- function()
{
   message(sprintf("--- test_splitExpressionMatrixByMutationStatusAtRSID_rosmap"))

     #-----------------------------------------------
     # use the eqtl-optimized matrix from sage
     #-----------------------------------------------

   variant.best <- "rs1136224"   # splits mtx.rna with t.test p.value < 1e-17
   variant.worst <- "rs352680"   # splits mtx.rna with t.test p.value > 0.99
   mtx.rna <- etx$get.rna.matrix("sage-eqtl-rosmap")

   targetGene <- "NDUFS2"
   etx <- EndophenotypeExplorer$new(targetGene, "hg38")

   x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx.rna, variant.best, study.name="rosmap")
   checkEquals(x$genotypes.vcf, list(wt=820, mut=324, het=296, hom=28))
   checkEquals(x$genotypes.rna, list(wt=388, mut=165, het=146, hom=19))

   ndufs2.wt <- as.numeric(x$wt["NDUFS2",])
   ndufs2.mut <- as.numeric(x$mut["NDUFS2",])
   checkTrue(t.test(ndufs2.wt, ndufs2.mut)$p.value < 1e-16)

   x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx.rna, variant.worst, study.name="rosmap")
   ndufs2.wt <- as.numeric(x$wt["NDUFS2",])
   ndufs2.mut <- as.numeric(x$mut["NDUFS2",])
   checkTrue(t.test(ndufs2.wt, ndufs2.mut)$p.value > 0.92)

   x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx.rna, "rs4575098", study.name="rosmap")
   ndufs2.wt <- as.numeric(x$wt["NDUFS2",])
   ndufs2.mut <- as.numeric(x$mut["NDUFS2",])
   checkTrue(t.test(ndufs2.wt, ndufs2.mut)$p.value < 6e-07)

   hsf2.wt <- as.numeric(x$wt["HSF2",])
   hsf2.mut <- as.numeric(x$mut["HSF2",])

   tfs.all <- get(load("~/github/MotifDb/inst/extdata/tfs/tfs-1683-lambert.RData"))$Gene
   tf.candidates <- intersect(tfs.all, rownames(mtx.rna))
   length(tf.candidates)

   library(trena)
   mtx.rna.fixed <- mtx.rna
   mtx.rna.fixed[is.na(mtx.rna.fixed)] <- 0
   solver <- EnsembleSolver(mtx.rna.fixed,
                           targetGene="NDUFS2",
                           candidateRegulators=tf.candidates,
                           solverNames=c("Spearman", "Pearson", "RandomForest", "xgboost"),
                           #solverNames=c("lasso", "Ridge", "Spearman", "Pearson", "RandomForest", "xgboost"),
			   geneCutoff=1.0)
   tbl.out <- run(solver)
   new.order <- order(abs(tbl.out$spearman), decreasing=TRUE)
   tbl.out <- tbl.out[new.order,]

   cor(hsf2.wt, ndufs2.wt, method="spearman", use="pairwise.complete")
   cor(hsf2.mut, ndufs2.mut, method="spearman", use="pairwise.complete")

} # test_splitExpressionMatrixByMutationStatusAtRSID_rosmap
#----------------------------------------------------------------------------------------------------
test_trenaScoreGenotypeStratifiedExpression <- function()
{
   message(sprintf("--- test_trenaScoreGenotypeStratifiedExpression"))

   dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
   file.cer <- "mtx.mayo.cer.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x255.RData"
   mtx.cer <- get(load(file.path(dir, file.cer)))

   tfs <- get(load(system.file(package="EndophenotypeExplorer", "extdata", "expression",
                               "tfs-44-for-PTK2B-mayo.RData")))

   checkEquals(length(tfs), 44)
   targetGene <- "PTK2B"
   etx <- EndophenotypeExplorer$new(targetGene, "hg38")
   tbl.eQTL <- etx$getEQTLsForGene()
   tbl.eQTL.sig <- subset(tbl.eQTL, pvalue < 0.0007 & study=="ampad-mayo")
   new.order <- order(tbl.eQTL.sig$pvalue, decreasing=FALSE)
   tbl.eQTL.sig <- tbl.eQTL.sig[new.order,]
   rsids <- tbl.eQTL.sig$rsid   #    "rs1594829"  "rs6557914"  "rs12544446"
   checkEquals(rsids, c("rs1594829", "rs6557914", "rs12544446"))

     # do just the first:
     # chrom                 chr8
     # hg19              26206077
     # hg38              26348561
     # rsid             rs1594829
     # pvalue         0.000129674
     # ensg       ENSG00000120899
     # genesymbol           PTK2B
     # study           ampad-mayo
     # tissue                 cer
     # assay              unknown

   rsid <- rsids[1]
   x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx.cer, rsid, study.name="mayo")
   expected <- c("wt", "mut", "het", "hom", "genotypes.rna", "genotypes.vcf")
   checkTrue(all(expected %in% names(x)))
   checkEquals(x$genotypes.vcf, list(wt=204, mut=145, het=127, hom=18))
   checkEquals(x$genotypes.rna, list(wt=155, mut=100, het=87, hom=13))
   checkEquals(ncol(x$wt), 155)
   checkEquals(ncol(x$mut), 100)
   checkEquals(ncol(x$het), 87)
   checkEquals(ncol(x$hom), 13)

      # with the rnaseq matrix now split by genotype, build and compare
      # trena models with two different pairs

   xx.mut <- etx$trenaScoreGenotypeStratifiedExpression(x$wt, x$mut, targetGene, tfs)
   xx.hom <- etx$trenaScoreGenotypeStratifiedExpression(x$wt, x$hom, targetGene, tfs)
   save(xx.mut, xx.hom, file="xx-mut-hom-tmp.RData")
   expected <- c("trena.1", "trena.2", "bicor.delta", "rf.delta", "spear.delta",
                 "lasso.delta")
   checkTrue(all(expected %in% names(xx.mut)))
   checkTrue(all(expected %in% names(xx.hom)))

   checkEquals(lapply(xx.mut, dim),
               list(trena.1=c(39, 8),
                    trena.2=c(34, 8),
                    bicor.delta=c(42, 5),
                    rf.delta=c(42, 5),
                    spear.delta=c(42, 5),
                    lasso.delta=c(42, 5)))

   checkEquals(lapply(xx.hom, dim),
               list(trena.1=c(39, 8),
                    trena.2=c(35, 8),
                    bicor.delta=c(43, 5),
                    rf.delta=c(43, 5),
                    spear.delta=c(43, 5),
                    lasso.delta=c(43, 5)))

      # combine these various delta tables for wt vs hom, choosing
      # just the top 2, sort by tf, see if any duplicats popup
      # this is but an adhoc way to summarize these data, adequate
      # only for this quick test.
      # see etx$summarizeResults for a more general and more satisfying
      # approach
   tbl.top.hom <- do.call(rbind, lapply(xx.hom[grep("delta", names(xx.hom))],
                                    function(tbl.delta) head(tbl.delta, 2)))
   rownames(tbl.top.hom) <- NULL
   tbl.top.hom <- tbl.top.hom[order(tbl.top.hom$tf),]
   checkEquals(length(grep("TBR1", tbl.top.hom$tf)), 2)

   tbl.top.mut <- do.call(rbind, lapply(xx.mut[grep("delta", names(xx.mut))],
                                    function(tbl.delta) head(tbl.delta, 2)))
   rownames(tbl.top.mut) <- NULL
   tbl.top.mut <- tbl.top.mut[order(tbl.top.mut$tf),]
   checkEquals(length(grep("TEAD1", tbl.top.mut$tf)), 3)

} # test_trenaScoreGenotypeStratifiedExpression
#----------------------------------------------------------------------------------------------------
test_summarizeStratifiedModels <- function()
{
   message(sprintf("--- test_summarizeStratifiedModels"))

   dir <- "~/github/TrenaProjectAD/prep/rna-seq-counts-from-synapse/eqtl"
   file.cer <- "mtx.mayo.cer.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x255.RData"
   mtx.cer <- get(load(file.path(dir, file.cer)))
   tfs <- get(load(system.file(package="EndophenotypeExplorer", "extdata", "expression",
                               "tfs-44-for-PTK2B-mayo.RData")))
   checkEquals(length(tfs), 44)
   targetGene <- "PTK2B"
   etx <- EndophenotypeExplorer$new(targetGene, "hg38")
   rsid <- "rs1594829"
   x <- etx$splitExpressionMatrixByMutationStatusAtRSID(mtx.cer, rsid, study.name="mayo")
   xx.mut <- etx$trenaScoreGenotypeStratifiedExpression(x$wt, x$mut, targetGene, tfs)
   xx.hom <- etx$trenaScoreGenotypeStratifiedExpression(x$wt, x$hom, targetGene, tfs)

   tbl.summary <- etx$summarizeStratifiedModels(xx.mut, "spearman")


} # test_summarizeStratifiedModels
#----------------------------------------------------------------------------------------------------
test_subsetAndRelabelGenoMatrixByPatientIDs <- function()
{
    message(sprintf("--- test_subsetAndRelabelGenoMatrixByPatientIDs"))

    etx <- EndophenotypeExplorer$new("NDUFS2", "hg38", initialize.snpLocs=TRUE)

    mtx.rna <- etx$get.rna.matrix("sage-eqtl-rosmap")
    variant.best  <- "rs1136224"   # from ampad eqtls
    variant.worst <- "rs352680"

    mtx.geno.1 <- etx$getGenoMatrixByRSID(variant.best)
    mtx.geno.pt.rosmap <- etx$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno.1, "rosmap")

    mtx.geno.pt.rosmap <- etx$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno.1, studyName="rosmap")

    checkEquals(dim(mtx.geno.1), c(1, 1894))
    checkEquals(dim(mtx.geno.pt.rosmap), c(1, 1151))

    tbl.map <- etx$getIdMap()

    expected.geno <-  c("0/0", "0/0", "0/1", "1/1", "0/0")
    vcf.sample.ids <- c("SM-CTELZ", "SM-CTEM1", "SM-CTEM2", "SM-CTEM3", "SM-CTEM4")
    patient.ids <- c("R5413523", "R9641047", "R9254359", "R4146432", "R6414128")
    checkEquals(as.character(mtx.geno.1[1, vcf.sample.ids]), expected.geno)
    checkEquals(as.character(mtx.geno.pt.rosmap[1, patient.ids]), expected.geno)

} # test_subsetAndRelabelGenoMatrixByPatientIDs
#----------------------------------------------------------------------------------------------------
# reproduce, more or less, these first and last eQTL scores from AMPAD
#  chrom      hg19      hg38       rsid       pvalue            ensg genesymbol        study tissue   assay
#   chr1 161184097 161214307  rs1136224 1.285003e-18 ENSG00000158864     NDUFS2 ampad-rosmap  dlpfc unknown
#   chr1 160678721 160708931   rs352680 0.9999419 ENSG00000158864     NDUFS2 ampad-rosmap  dlpfc unknown
demo_variantSeparated.rna.t.test <- function()
{
    message(sprintf("--- test_subsetAndRelabelGenoMatrixByPatientIDs"))

    targetGene <- "NDUFS2"
    etx <- EndophenotypeExplorer$new(targetGene, "hg38", initialize.snpLocs=TRUE)

    mtx.rna <- etx$get.rna.matrix("sage-eqtl-rosmap")
    variant.best  <- "rs1136224"   # from ampad eqtls
    variant.worst <- "rs352680"

    mtx.geno.1 <- etx$getGenoMatrixByRSID(variant.best)
    mtx.geno.pt.rosmap <- etx$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno.1, "rosmap")

    wt.patients <- which(mtx.geno.pt.rosmap[1,] %in% "0/0")
    wt.patient.names <- colnames(mtx.geno.pt.rosmap)[wt.patients]
    checkEquals(length(wt.patient.names), 827)
    wt.patient.names <- intersect(colnames(mtx.rna), wt.patient.names)
    checkEquals(length(wt.patient.names), 388)

    mut.patients <- which(mtx.geno.pt.rosmap[1,] %in% c("0/1", "1/1"))
    mut.patient.names <- colnames(mtx.geno.pt.rosmap)[mut.patients]
    checkEquals(length(mut.patient.names), 324)
    mut.patient.names <- intersect(colnames(mtx.rna), mut.patient.names)
    checkEquals(length(mut.patient.names), 165)

    wt.vec <- mtx.rna[targetGene, wt.patient.names]
    mut.vec <- mtx.rna[targetGene, mut.patient.names]
    checkEqualsNumeric(t.test(wt.vec, mut.vec)$p.value, 5.31e-17, tol=1e15)

       #------------------------------------------------
       # does mutation status correspond to diagnosis?
       # from mayo/sage ampad eqtl 2020 nature paper
       #  Prior to RNA-seq normalization, we harmonized the LOAD
       #  definition across AMP-AD studies. AD controls were defined as
       #  patients with a low burden of plaques and tangles, as well as
       #  lack of evidence of cognitive impairment. For the ROSMAP
       #  study, we defined AD cases to be individuals with a Braak
       #  greater than or equal to 4, CERAD score less than or equal
       #  to 2, and a cognitive diagnosis of probable AD with no other
       #  causes (cogdx = 4), and controls to be individuals with Braak
       #  less than or equal to 3, CERAD score greater than or equal to
       #  3, and cognitive diagnosis of ‘no cognitive impairment’
       #  (cogdx = 1). For the Mayo Clinic study, we defined disease
       #  status based on neuropathology, where individuals with Braak
       #  score greater than or equal to 4 were defined to be AD cases,
       #  and individuals with Braak less than or equal to 3 were
       #  defined to be controls. Individuals not meeting “AD case” or
       #  “control” criteria were retained for analysis, and were
       #  categorized as “other” for the purposes of RNA-seq
       #  adjustment.

    tbl.pt <- etx$get.rosmap.patient.table(NA)
    checkEquals(dim(tbl.pt), c(3583, 18))
        # table(tbl.pt$braaksc)
        #  0   1   2   3   4   5   6
        #  17 104 152 383 510 408  23
        # table(tbl.pt$ceradsc)
        #   1   2   3   4
        # 532 560 141 364
        # table(tbl.pt$cogdx)
        #   1   2   3   4   5   6
        # 586 404  33 674  94  30

     tbl.pt.ad <- subset(tbl.pt, braaksc >=4 & ceradsc <= 2 & cogdx==4)
     checkEquals(dim(tbl.pt.ad), c(445, 18))
     pts.ad <- unique(tbl.pt.ad$individualID)
     pts.ad <- intersect(pts.ad, colnames(mtx.geno.pt.rosmap))
     length(pts.ad)  # 312

     tbl.pt.ctl <- subset(tbl.pt, braaksc <=3 & ceradsc >= 3 & cogdx==1)
     checkEquals(dim(tbl.pt.ctl), c(203, 18))
     pts.ctl <- unique(tbl.pt.ctl$individualID)
     pts.ctl <- intersect(pts.ctl, colnames(mtx.geno.pt.rosmap))
     length(pts.ctl)  # 150

        # no phenotype difference
     round(as.numeric(table(mtx.geno.pt.rosmap[, pts.ctl]))/150, digits=2) # [1] 0.69 0.30 0.01
     round(as.numeric(table(mtx.geno.pt.rosmap[, pts.ad]))/312, digits=2)  # [1] 0.71 0.27 0.02

} # demo_variantSeparated.rna.t.test
#----------------------------------------------------------------------------------------------------
exploreGenotypeExpressionPhenotype <- function()
{
   #rsid <-

} # explore genotype/expression/AD-phenotype
#----------------------------------------------------------------------------------------------------
if(!interactive()){
   init.snpLocs()
   runTests()
   }
