# EndophenotypeExplorer
#----------------------------------------------------------------------------------------------------
#' @importFrom RPostgreSQL dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @import RPostgreSQL
#' @import org.Hs.eg.db
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import R6
#' @import VariantAnnotation
#' @import GenomicRanges
#' @importFrom rtracklayer liftOver import.chain
#' @import httr
#' @import jsonlite
#'
#' @title EndophenotypeExplorer
#------------------------------------------------------------------------------------------------------------------------
#' @name EndophenotypeExplorer-class
#' @rdname EndophenotypeExplorer-class
#' @aliases EndophenotypeExplorer
#'
# library(R6)
# library(VariantAnnotation)
# library(RPostgreSQL)
#----------------------------------------------------------------------------------------------------
#' R6 Class for exploring associations between eQTLs, variants, gene expression and gene regulation
#'
#' An endophenotypeExplore has access to eQTLs, variant calls, and more
#'
#' @export
#'
#----------------------------------------------------------------------------------------------------
EndophenotypeExplorer = R6Class("EndophenotypeExplorer",

    #--------------------------------------------------------------------------------
    private = list(target.gene=NULL,
                   chromosome=NULL,
                   vcf.url=NULL,
                   geneReg.db=NULL,
                   default.genome=NULL,
                   tbl.clinical=NULL,
                   tbl.idMap=NULL,
                   tbl.biospecimen.rosmap=NULL,
                   tbl.biospecimen.mayo=NULL,
                   tbl.biospecimen.sinai=NULL,
                   tbl.clinical.rosmap=NULL,
                   tbl.clinical.sinai=NULL,
                   tbl.clinical.mayo=NULL,
                   tbl.gwas.38=NULL,
                   tbl.gwas.38.associated=NULL,
                   verbose=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(

      #' @description
      #' Create a new EndophenotypeExplorer object
      #' @param target.gene  Gene of interest.
      #' @param default.genome UCSC code, either `hg19` or `hg38`.
      #' @param vcf.url https endpoint from serving indexed vcf files
      #' @return A new `EndophenotypeExplorer` object.

        initialize = function(target.gene, default.genome, verbose=FALSE){
            private$target.gene <- target.gene
            private$default.genome <- default.genome
            private$chromosome <- self$identifyTargetGeneChromosome(target.gene)
            private$vcf.url <- self$setupVcfURL(private$chromosome)
            private$verbose <- verbose
            self$setupClinicalData()
            self$setupGWASData()
            },

        setupVcfURL = function(chromosome){
           vcf.base.url <- "https://igv-data.systemsbiology.net/static/ampad"
           vcf.directory <- "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_recalibrated_variants"
           sprintf("%s/%s/%s.vcf.gz", vcf.base.url, vcf.directory, chromosome)
           },

        getVcfUrl = function(){
           private$vcf.url
           },

        getGenoMatrix = function(chrom, start, end){
            roi <- GRanges(seqnames=chrom, IRanges(start=start, end=end))
            x <- readVcf(private$vcf.url, private$default.genome, roi)
            stopifnot("GT" %in% names(geno(x)))
            invisible(geno(x)$GT)
            },

        getIdMap = function(){
           private$tbl.idMap
           },

        mapSampleIdToPatientAndCohort = function(sampleID){
           return(subset(private$tbl.idMap, vcf==sampleID))
           },

        setupClinicalData = function(){
               # first the per-sample biospecimen data
            dir <- system.file(package="EndophenotypeExplorer", "extdata", "biospecimen")
            private$tbl.biospecimen.rosmap <- read.table(file.path(dir, "ROSMAP_biospecimen_metadata.csv"),
                                                         sep=",", header=TRUE, as.is=TRUE)
            private$tbl.biospecimen.sinai <- read.table(file.path(dir, "MSBB_biospecimen_metadata.csv"),
                                                        sep=",", header=TRUE, as.is=TRUE)
            private$tbl.biospecimen.mayo <- read.table(file.path(dir, "MayoRNAseq_biospecimen_metadata.csv"),
                                                       sep=",", header=TRUE, as.is=TRUE)
              # now the per-patient clinical data
            dir <- system.file(package="EndophenotypeExplorer", "extdata", "clinical")
            private$tbl.clinical.rosmap <- read.table(file.path(dir, "ROSMAP_clinical.csv"),
                                                      sep=",", header=TRUE, as.is=TRUE)
            private$tbl.clinical.sinai <- read.table(file.path(dir, "MSBB_individual_metadata.csv"),
                                                      sep=",", header=TRUE, as.is=TRUE)
            private$tbl.clinical.mayo <- read.table(file.path(dir, "MayoRNAseq_individual_metadata.csv"),
                                                      sep=",", header=TRUE, as.is=TRUE)
            dir <- system.file(package="EndophenotypeExplorer", "extdata", "idMapping")
               # the crucial sample-to-patient mapping
               # this was difficult to obtain.  see the mapToPatientID function in
               # ~/github/TrenaProjectAD/explore/ampad.eQTLS/ldlr.R
               # todo: move this code to a prep directory in this package
            if(private$verbose) message(sprintf("--- about to load id mapping file"))
            if(private$verbose) message(sprintf("    dir: %s", dir))
            full.path <- file.path(dir, "tbl.vcfToPatientIDs.RData")
            if(private$verbose) message(sprintf("    full.path: %s", full.path))
            private$tbl.idMap <- get(load(file.path(dir, "tbl.vcfToPatientIDs.RData")))
            },

        getClinicalTable = function(){
            private$tbl.clinical
            },

        setupGWASData = function(){
           dir <- system.file(package="EndophenotypeExplorer", "extdata", "gwas")
           full.path <- file.path(dir, "tbl.posthuma-38-loci-curated.RData")
           private$tbl.gwas.38 <- get(load(full.path))
           full.path <- file.path(dir, "tbl.posthuma-38-geneAssociations-curated-3828x12.RData")
           private$tbl.gwas.38.associated <- get(load(full.path))
           },

        getGWASTables = function(){
           list(gwas.38 = private$tbl.gwas.38,
                gwas.38.assoc = private$tbl.gwas.38.associated)
           },

        identifyTargetGeneChromosome = function(target.gene){
           suppressMessages({
             entrezID <- select(org.Hs.eg.db, keys=target.gene, keytype="SYMBOL", columns="ENTREZID")$ENTREZID
             stopifnot(nchar(entrezID) > 1)
             select(TxDb.Hsapiens.UCSC.hg38.knownGene, keys=entrezID, keytype="GENEID",
                    columns="TXCHROM")$TXCHROM
             })
           },

        get.rosmap.patient.table = function(patientID){
            if(is.na(patientID))
               return(private$tbl.clinical.rosmap)
            subset(private$tbl.clinical.rosmap, individualID==patientID)
            },

        get.sinai.patient.table = function(patientID){
            if(is.na(patientID))
               return(private$tbl.clinical.sinai)
            subset(private$tbl.clinical.sinai, individualID==patientID)
            },

        get.mayo.patient.table = function(patientID){
            if(is.na(patientID))
               return(private$tbl.clinical.mayo)
            subset(private$tbl.clinical.mayo, individualID==patientID)
            },

        vcfSampleID.to.clinicalTable = function(sampleID){
            if(!sampleID %in% private$tbl.idMap$vcf)
                return(data.frame())
            tbl <- subset(private$tbl.idMap, vcf==sampleID)
            patientID <- tbl$patient
            cohort <- tbl$cohort
            tbl.patient <- switch(cohort,
              "mayo"   = self$standardizeMayoPatientTable(subset(private$tbl.clinical.mayo,
                                                          individualID==patientID)),
              "rosmap" = self$standardizeRosmapPatientTable(subset(private$tbl.clinical.rosmap,
                                                            individualID==patientID)),
              "sinai"  = self$standardizeSinaiPatientTable(subset(private$tbl.clinical.sinai,
                                                           individualID==patientID))
              )
            tbl.patient$cohort <- cohort
            tbl.patient$sampleID <- sampleID
            coi <- c("patientID", "sampleID", "cohort", "study","sex","ethnicity","apoeGenotype",
                     "braak","cerad","cogdx","pmi", "ageAtDeath")
            tbl.patient[, coi]
            },

        standardizeMayoPatientTable = function(tbl){
           tbl$ageDeath <- as.numeric(sub("+", "", tbl$ageDeath, fixed=TRUE))
           coi <- c("individualID", "individualIdSource", "sex", "ethnicity", "apoeGenotype",
                    "Braak", "CERAD", "pmi", "ageDeath")
           tbl.0 <- tbl[, coi]
           tbl.0$cogdx <- NA
           standard.names <- c("patientID", "study", "sex","ethnicity","apoeGenotype","braak","cerad",
                                "pmi", "ageAtDeath", "cogdx")
           colnames(tbl.0) <- standard.names
           tbl.0
           },

        standardizeRosmapPatientTable = function(tbl){
           tbl$age_death <- as.numeric(sub("+", "", tbl$age_death, fixed=TRUE))
           coi <- c("individualID", "Study", "msex", "race", "apoe_genotype", "braaksc", "ceradsc",
                    "cogdx", "pmi", "age_death")
           tbl.0 <- tbl[, coi]
           standard.names <- c("patientID", "study", "sex","ethnicity","apoeGenotype","braak","cerad",
                                "cogdx", "pmi", "ageAtDeath")
           colnames(tbl.0) <- standard.names
           tbl.0
           },

        standardizeSinaiPatientTable = function(tbl){
           tbl$ageDeath <- as.numeric(sub("+", "", tbl$ageDeath, fixed=TRUE))
           coi <- c("individualID","individualIdSource","sex","ethnicity","apoeGenotype",
                    "Braak","CERAD","pmi","ageDeath")
           tbl.0 <- tbl[, coi]
           tbl.0$cogdx <- NA
           standard.names <- c("patientID", "study", "sex","ethnicity","apoeGenotype","braak","cerad",
                                "pmi", "ageAtDeath", "cogdx")
           colnames(tbl.0) <- standard.names
           tbl.0
           },

        getAggregatedAlleleFrequencies = function(rsid){

            uri <- sprintf("https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/%s/frequency",
                  sub("rs", "", rsid))
            response <- GET(uri)
            suppressMessages(x <- fromJSON(httr::content(response, as="text"))$results)
            ref <- x[[1]]$ref
            counts.list <- x[[1]]$counts

            population.codes <- names(counts.list[[1]][[1]])
            tbls <- list()
            for(pop in population.codes){
                pop.record <- counts.list[[1]][[1]][[pop]]
                tbl <- data.frame(population=pop, ref=ref, rsid=rsid, stringsAsFactors=FALSE)
                alleles <- names(pop.record)
                counts <- as.integer(pop.record)
                tbl[1, alleles] <- counts
                tbls[[pop]] <- tbl
            }
            tbl.alfa <- do.call(rbind, tbls)
            rownames(tbl.alfa) <- NULL

            f <- system.file(package="EndophenotypeExplorer", "extdata", "gwas", "populationCodes.tsv")
            stopifnot(file.exists(f))
            tbl.popCodes <- read.table(f, header=TRUE, as.is=TRUE, sep="\t")

            tbl.alfa <- merge(tbl.alfa, tbl.popCodes, by.x="population", by.y="biosampleID")
            allele.colnames <- setdiff(colnames(tbl.alfa), c("population", "ref", "name", "rsid"))
            tbl.alfa <- tbl.alfa[, c("rsid", "ref", "name", allele.colnames)]
            deleter <- grep("population", colnames(tbl.alfa))
            if(length(deleter) == 1)
                tbl.alfa <- tbl.alfa[, -(deleter)]
            coi <- c("rsid", "ref", "name", allele.colnames)
            tbl.alfa <- tbl.alfa[, coi]
            colnames(tbl.alfa)[grep("name", colnames(tbl.alfa))] <- "population"
            tbl.alfa$total <- rowSums(tbl.alfa[, allele.colnames])

            for(allele in allele.colnames){
                allele.freq <- 100 * tbl.alfa[[allele]]/tbl.alfa$total
                new.colname <- sprintf("%s.freq", allele)
                tbl.alfa[[new.colname]] <- allele.freq
            } # for allele
           tbl.alfa
           } # getAggregatedAlleleFrequencies

       ) # public

    ) # class EndophenotypeExplorer

#--------------------------------------------------------------------------------


