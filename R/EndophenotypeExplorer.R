# EndophenotypeExplorer
#----------------------------------------------------------------------------------------------------
#' @importFrom RPostgreSQL dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @import RPostgreSQL
#' @import org.Hs.eg.db
#' @importFrom BSgenome snpsById snpsByOverlaps
#' @import SNPlocs.Hsapiens.dbSNP144.GRCh37
#' @import SNPlocs.Hsapiens.dbSNP151.GRCh38
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
           vcf.directory <- "NIA-1898"
           sprintf("%s/%s/%s.vcf.gz", vcf.base.url, vcf.directory, chromosome)
           },

        getVcfUrl = function(){
           private$vcf.url
           },

        getGenoMatrix = function(chrom, start, end){
            printf("entering getGenoMatrix with vectors")
            chromosomesInProperFormat <- !grepl("chr", chrom[1])
            stopifnot(chromosomesInProperFormat)
            roi <- GRanges(seqnames=chrom, IRanges(start=start, end=end))
            x <- readVcf(private$vcf.url, private$default.genome, roi)
            stopifnot("GT" %in% names(geno(x)))
            mtx <- geno(x)$GT
            invisible(mtx)
            },

        getGenoMatrixByRSID = function(rsids){
            tbl.locs <- self$rsidToLoc(rsids)
            mtx <- self$getGenoMatrix(tbl.locs$chrom, tbl.locs$hg19, tbl.locs$hg19)
            #roi <- with(tbl.locs, GRanges(seqnames=chrom, IRanges(start=hg19, end=hg19)))
            #x <- readVcf(private$vcf.url, private$default.genome, roi)
            #stopifnot("GT" %in% names(geno(x)))
            #mtx <- geno(x)$GT
            invisible(mtx)
            },

        getIdMap = function(){
           private$tbl.idMap
           },

        rsidToLoc = function(rsids){
            rsids <- grep("^rs", rsids, value=TRUE)
            gr.hg19 <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, rsids, ifnotfound="drop")
            tbl.hg19 <- as.data.frame(gr.hg19)
            colnames(tbl.hg19)[1] <- "chrom"
            colnames(tbl.hg19)[2] <- "hg19"
            colnames(tbl.hg19)[4] <- "rsid"
            tbl.hg19$chrom <- as.character(tbl.hg19$chrom)
            gr.hg38 <- snpsById(SNPlocs.Hsapiens.dbSNP151.GRCh38, rsids, ifnotfound="drop")
            tbl.hg38 <- as.data.frame(gr.hg38)
            colnames(tbl.hg38)[2] <- "hg38"
            colnames(tbl.hg38)[4] <- "rsid"
            tbl.out <- merge(tbl.hg19[, c("chrom", "hg19", "rsid")], tbl.hg38[, c("hg38", "rsid")], by="rsid")
            tbl.out[, c("chrom", "hg19", "hg38", "rsid")]
            },

        locsToRSID = function(locs, genome){
            chroms <- unlist(lapply(strsplit(locs, ":"), "[", 1))
            loc.strings <- unlist(lapply(strsplit(locs, ":"), "[", 2))
            locs <- as.integer(sub("_.*$", "", loc.strings))
            snplocs <- switch(genome,
                   "hg19" = SNPlocs.Hsapiens.dbSNP144.GRCh37,
                   "hg38" = SNPlocs.Hsapiens.dbSNP151.GRCh38
                   )
            gr <- GRanges(seqnames=chroms, IRanges(start=locs, end=locs))
            gr.snps <- snpsByOverlaps(snplocs, gr)
            tbl.rsids <- as.data.frame(gr.snps)[, c("seqnames", "pos", "RefSNP_id")]
            colnames(tbl.rsids)[1] <- "chrom"
            colnames(tbl.rsids)[2] <- "loc"
            colnames(tbl.rsids)[3] <- "rsid"
            tbl.rsids$chrom <- as.character(tbl.rsids$chrom)
            tbl.rsids$signature <- sprintf("%s:%s", tbl.rsids$chrom, tbl.rsids$loc)
            sig <- sprintf("%s:%s", chroms, locs)
            tbl.all <- data.frame(chrom=chroms, loc=locs,
                                  signature=sig,
                                  stringsAsFactors=FALSE)
            tbl.new <- merge(tbl.all, tbl.rsids[, c("rsid", "signature")], by="signature", all.x=TRUE)
            failures <- which(is.na(tbl.new$rsid))
            length(failures)
            tbl.new$rsid[failures] <- tbl.new$signature[failures]
            tbl.new$rsid
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
            min.freqs <- unlist(lapply(seq_len(nrow(tbl.alfa)), function(r){
                freqs <- tbl.alfa[r,grep("freq", colnames(tbl.alfa))]
                freqs <- freqs[freqs>0]
                min(freqs, na.rm=TRUE)
                }))
            tbl.alfa$min.freq <- min.freqs
            tbl.alfa
            }, # getAggregatedAlleleFrequencies

        getTissueExpression = function(geneSymbol, tissue){
           stopifnot(tolower(tissue)=="brain")
           f <- system.file(package="EndophenotypeExplorer", "extdata", "gtex", "tissues.RData")
           supportedTissues <- get(load(f))
           f <- system.file(package="EndophenotypeExplorer", "extdata", "gtex", "tbl.gencode-v26-GRCh38.RData")
           tbl.gencode <- get(load(f))
           stopifnot(geneSymbol %in% tbl.gencode$gene_name)
           gtex.tissues <- grep("brain", supportedTissues, ignore.case=TRUE, value=TRUE)
           gencode.id <- subset(tbl.gencode, gene_name==geneSymbol)$gene_id[1]
           tissues.string <- paste(gtex.tissues, collapse=",")
           url.0 <- "https://gtexportal.org/rest/v1/expression/medianGeneExpression"
           url.1 <- "datasetId=gtex_v8"
           url.2 <- sprintf("gencodeId=%s", gencode.id)
           url.3 <- sprintf("tissueSiteDetailId=%s", tissues.string)
           url.4 <- "format=json"
           url <- sprintf("%s?%s&%s&%s&%s", url.0, url.1, url.2, url.3, url.4)
           x <- content(GET(url))
           tbls <- lapply(x[[1]], as.data.frame)
           tbl <- do.call(rbind, tbls)
           tbl
           }

       ) # public

    ) # class EndophenotypeExplorer

#--------------------------------------------------------------------------------


