# EndophenotypeExplorer
#----------------------------------------------------------------------------------------------------
#' @importFrom RPostgreSQL dbConnect dbListTables dbGetQuery dbListConnections dbDisconnect
#' @import RPostgreSQL
#' @import org.Hs.eg.db
#' @importFrom BSgenome snpsById snpsByOverlaps
#' @import SNPlocs.Hsapiens.dbSNP144.GRCh37
#' @import SNPlocs.Hsapiens.dbSNP155.GRCh38
#' @import biomaRt
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import R6
#' @import VariantAnnotation
#' @import GenomicRanges
#' @importFrom rtracklayer liftOver import.chain
#' @import httr
#' @import jsonlite
#' @import trena
#'
#' @title EndophenotypeExplorer
#------------------------------------------------------------------------------------------------------------------------
#' @name EndophenotypeExplorer-class
#' @rdname EndophenotypeExplorer-class
#' @aliases EndophenotypeExplorer
#----------------------------------------------------------------------------------------------------
#'@description
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
                   vcf.project=NULL,
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
                   tbl.eqtls=NULL,
                   standard.clinical.columns=NULL,
                   snpMart.hg38=NULL,
                   snpMart.hg19=NULL,
                   expressionMatrixCodes=NULL,
                   verbose=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(

      #' @description
      #' Create a new EndophenotypeExplorer object
      #' @param target.gene  Gene of interest.
      #' @param default.genome UCSC code, either `hg19` or `hg38`.
      #' @param vcf.project character, either "ADNI" or "AMPAD"
      #' @param chromosome default NA, may be unusual cases where we supply it (C2orf90/ECRG4) for instance
      #' @param verbose logical
      #' @param initialize.snpLocs logical, force this ~60 second process at startup
      #' @parma defer.setupClinicalData.toSupportTesting logical, default FALSE
      #' @return A new `EndophenotypeExplorer` object.

        initialize = function(target.gene, default.genome,
                              vcf.project,
                              chromosome=NA,
                              verbose=FALSE,
                              initialize.snpLocs=FALSE,
                              defer.setupClinicalData.toSupportTesting=FALSE){
            stopifnot(vcf.project %in% c("AMPAD", "ADNI"))
            if(!default.genome %in% c("hg19", "hg38")){
               message(sprintf("error: genome '%s' not supported", default.genome))
               stop()
               }
            private$vcf.project <- vcf.project
            private$default.genome = default.genome
            self$setTargetGene(target.gene, default.genome, chromosome)
            private$verbose <- verbose
            private$standard.clinical.columns <- c("patientID", "study", "sex","ethnicity",
                                                   "apoeGenotype","braak","cerad", "pmi",
                                                   "ageAtDeath", "cogdx")
            self$setupGWASData()
            if(initialize.snpLocs){
                message(sprintf("initializing hg19 and hg38 snpLocs, may take a minute"))
                t1 <- system.time(x <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, "rs769450"))
                t2 <- system.time(x <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, "rs769450"))
                message(sprintf("dbSNP144 hg19: %5.2f", t2[["elapsed"]]))
                message(sprintf("dbSNP155 hg38: %5.2f", t1[["elapsed"]]))
                }
            if(!defer.setupClinicalData.toSupportTesting)
                self$setupClinicalData()
            private$expressionMatrixCodes <-
                list(
                 "max-tcx"="mayo.tcx.robinson.normalized.PMI-age-cellType-covariate-collected.15201x262.RData",
                 "sage-eqtl-cer"="mtx.mayo.cer.eqtl-optimized-geneSymbols-sampleIDs-17009x261.RData",
                 "sage-eqtl-tcx"="mtx.mayo.tcx.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x257.RData",
                 "old-mayo-tcx"="mayo.tcx.16969x262.covariateCorrection.log+scale.RData",
                 "old-mayo-cer"="cerebellum.15167x263.RData",
                 "old-rosmap"="rosmap.14235x632.RData",
                 "sage-eqtl-rosmap"="mtx.rosmap.rnaseq-residual-eqtl-geneSymbols-patients-15582x632.RData",
                 "sage-counts-rosmap"="mtx.rosmap.rnaseq-counts-geneSymbols-patients-15582x632.RData",
                 "max-rosmap"="ROSMAP_rnaseq-ncsNormalized-covCorrectedForStudyPMIandAgeAtDeath-10oct21-21667x631.RData",
                 "tcx-unknown"="temporalCortex.15167x264.RData",
                 "GTEx_V8.Brain_Amygdala"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Amygdala.RData",
                 "GTEx_V8.Brain_Anterior_cingulate_cortex_BA24"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Anterior_cingulate_cortex_BA24.RData",
                 "GTEx_V8.Brain_Caudate_basal_ganglia"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Caudate_basal_ganglia.RData",
                 "GTEx_V8.Brain_Cerebellar_Hemisphere"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Cerebellar_Hemisphere.RData",
                 "GTEx_V8.Brain_Cerebellum"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Cerebellum.RData",
                 "GTEx_V8.Brain_Cortex"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Cortex.RData",
                 "GTEx_V8.Brain_Frontal_Cortex_BA9"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Frontal_Cortex_BA9.RData",
                 "GTEx_V8.Brain_Hippocampus"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Hippocampus.RData",
                 "GTEx_V8.Brain_Hypothalamus"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Hypothalamus.RData",
                 "GTEx_V8.Brain_Nucleus_accumbens_basal_ganglia"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Nucleus_accumbens_basal_ganglia.RData",
                 "GTEx_V8.Brain_Putamen_basal_ganglia"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Putamen_basal_ganglia.RData",
                 "GTEx_V8.Brain_Spinal_cord_cervical_c-1"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Spinal_cord_cervical_c-1.RData",
                 "GTEx_V8.Brain_Substantia_nigra"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Brain_Substantia_nigra.RData",
                 "GTEx_V8.Whole_Blood"="~/github/TrenaProjectAD/inst/extdata/expression/gtex.v8/Whole_Blood.RData"
                 )
            },

         #' @description
         #' the targetGene and genome are widely referenced - this member function allows you
         #' to change your initial choices
         #' @param target.gene  Gene of interest.
         #' @param genome UCSC code, either `hg19` or `hg38`.
         #' @param chromosome default NA, causing a lookup
        setTargetGene = function(targetGene, genome, chromosome=NA){
            private$target.gene <- targetGene
            private$default.genome <- genome
            if(is.na(targetGene)){
              private$chromosome <- NA
              private$vcf.url <- NA
            } else {
               if(!is.na(chromosome)){
                   private$chromosome <- chromosome
                   } else {
                      private$chromosome <- self$identifyTargetGeneChromosome(targetGene)
                      }
               private$vcf.url <- self$setupVcfURL(private$chromosome)
               }
            },

         #' @description
         #' the target chromosome and genome are widely referenced - this member function allows you
         #' to change your initial choices
         #' @param target.chromosome
         #' @param genome UCSC code, either `hg19` or `hg38`.
        setTargetChromosome = function(target.chromosome, genome){
            private$default.genome <- genome
            private$chromosome <- target.chromosome
            private$vcf.url <- self$setupVcfURL(private$chromosome)
            },

        get.rna.matrix.codes=function(){
           return(private$expressionMatrixCodes)
           },

        get.rna.matrix=function(code) {
            stopifnot(Sys.info()[["user"]] %in% c("paul", "pshannon"))  # only works for me, for now
            mtx.choices <- names(private$expressionMatrixCodes)
            stopifnot(code %in% mtx.choices)

            if(grepl("^GTEx_V8", code)){
               path <- private$expressionMatrixCodes[[code]]
               mtx.rna <- get(load(path))
               }

            if(code=="max-rosmap"){
               data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression/maxNormalizations"
               filename <- "ROSMAP_rnaseq-ncsNormalized-covCorrectedForStudyPMIandAgeAtDeath-10oct21-21667x631.RData"
               mtx.rna <- get(load(file.path(data.dir, filename)))
               }

            if(code=="max-tcx"){
               data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression/maxNormalizations"
               filename <- "mayo.tcx.robinson.normalized.PMI-age-cellType-covariate-collected.15201x262.RData"
               mtx.rna <- get(load(file.path(data.dir, filename)))
               }

            if(code=="sage-eqtl-cer"){
               data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression/sage.eqtl.optimized"
               filename <- "mtx.mayo.cer.eqtl-optimized-geneSymbols-sampleIDs-17009x261.RData"
               mtx.rna <- get(load(file.path(data.dir, filename)))
               }

            if(code=="sage-eqtl-tcx"){
                data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression/sage.eqtl.optimized"
                filename <- "mtx.mayo.tcx.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x257.RData"
                mtx.rna <- get(load(file.path(data.dir, filename)))
            }

            if(code=="old-mayo-tcx"){
                data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression"
                filename <- "mayo.tcx.16969x262.covariateCorrection.log+scale.RData"
                mtx.rna <- get(load(file.path(data.dir, filename)))
            }

            if(code=="old-mayo-cer"){
                data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression"
                filename <- "cerebellum.15167x263.RData"
                mtx.rna <- get(load(file.path(data.dir, filename)))
            }

            if(code=="old-rosmap"){
                data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression"
                filename <- "rosmap.14235x632.RData"
                mtx.rna <- get(load(file.path(data.dir, filename)))
            }

            if(code=="sage-eqtl-rosmap"){
                data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression/sage.eqtl.optimized"
                filename <- "mtx.rosmap.rnaseq-residual-eqtl-geneSymbols-patients-15582x632.RData"
                mtx.rna <- get(load(file.path(data.dir, filename)))
            }
            if(code=="sage-counts-rosmap"){
                data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression/sage.eqtl.optimized"
                filename <- "mtx.rosmap.rnaseq-counts-geneSymbols-patients-15582x632.RData"
                mtx.rna <- get(load(file.path(data.dir, filename)))
            }

            if(code=="tcx-unknown"){
                data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression"
                filename <- "temporalCortex.15167x264.RData"
                mtx.rna <- get(load(file.path(data.dir, filename)))
            }

            invisible(mtx.rna)
            }, # get.rna.matrix
#----------------------------------------------------------------------------------------------------

         #' @description
         #' our genotypes are stored per-chromosome in vcf files.  this
         #' function uses the current chromosome of interest, to provide
         #' otherwise simple access to genotypes of the current gene and its
         #' putative regulatory regions
         #' @param chromosome
         #' @returns the url
        setupVcfURL = function(chromosome){
           url <- "failure"
           if(private$vcf.project == "AMPAD"){
              vcf.base.url <- "https://igv-data.systemsbiology.net/ampad"
              vcf.directory <- "NIA-1898"
              url <- sprintf("%s/%s/%s.vcf.gz", vcf.base.url, vcf.directory, chromosome)
              }
           if(private$vcf.project == "ADNI"){
              vcf.base.url <- "https://igv-data.systemsbiology.net/nfs/adni"
              if(chromosome == "M"){
                 vcf.filename <- "adni_mito_genomes_20170201.vcf.gz"
                 }
              if(chromosome != "M"){
                 template <- "gcad.qc.wgs.%s.4789.GATK.2018.09.17.v2.biallelic.genotypes.ALL.vcf.gz"
                 vcf.filename <- sprintf(template, chromosome)
                 }
              url <- sprintf("%s/%s", vcf.base.url, vcf.filename)
              } # ADNI
           return(url)
           },

        #' @description
        #' easy access to the current value
        #' returns url in character form
        getVcfUrl = function(){
           private$vcf.url
           },


        #' @description
        #' genotypes for all samples in the specified locus
        #' @param chrom character
        #' @param start numeric
        #' @param end numeric
        #' returns matrix, location by sample
        getGenoMatrix = function(chrom, start, end){
            fixChromFormat <- function(chrom){
              stopifnot(private$vcf.project %in% c("ADNI", "AMPAD"))
              stopifnot(all(chrom == chrom[1]))
              chrom.fixed  <- chrom[1]
              if(private$vcf.project == "ADNI"){
                  if(!grepl("^chr", chrom[1]))
                     chrom.fixed <- paste0("chr", chrom[1])
                 } # ADNI
              if(private$vcf.project == "AMPAD"){
                  if(grepl("^chr", chrom[1]))
                     chrom.fixed <- sub("chr", "", chrom[1])
                 } # AMPAD
               return(rep(chrom.fixed, length(chrom)))
               } # fixChromFormat
            chrom <- fixChromFormat(chrom)
            roi <- sort(GRanges(seqnames=chrom, IRanges(start=start, end=end)))
            x <- readVcf(private$vcf.url, private$default.genome, roi)
            stopifnot("GT" %in% names(geno(x)))
            mtx <- geno(x)$GT
            invisible(mtx)
            },

        #' @description
        #' genotypes for all samples at the specified rsids
        #' @param rsids character
        #' returns matrix, location by sample

        getGenoMatrixByRSID = function(rsids){
            stopifnot(all(grepl("^rs", rsids)))
            message(sprintf("--- getGenoMatrixByRSID looking up locs: %s", Sys.time()))
            tbl.locs <- self$rsidToLoc(rsids)
            message(sprintf("--- getGenoMatrixByRSID looking up locs, done: %s", Sys.time()))
            if(nrow(tbl.locs) == 0){
                message(sprintf("failed to find chrom loc for rsids: %s",
                                paste(rsids, collapse=",")))
                return(NA)
                }
            start <- tbl.locs[, private$default.genome]
            end <- start
            mtx <- self$getGenoMatrix(tbl.locs$chrom, start, end)
            message(sprintf("--- getGenoMatrixByRSID, back from getGenoMatrix (locs): %s", Sys.time()))
            invisible(mtx)
            },

        get.ampad.EQTLsForGene = function(){
            if(grepl("hagfish", Sys.info()[["nodename"]])){
               suppressWarnings(db.access.test <- try(system("/sbin/ping -c 1 khaleesi", intern=TRUE, ignore.stderr=TRUE)))
               if(length(db.access.test) == 0){
                  message(sprintf(("khaleesi unreachable, no eQTLS available")))
                  return(data.frame())
                  } # length 0
               } # on macos hagfish
            db <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="genereg2021", host="khaleesi")
            query.string <- sprintf("select * from eqtls where genesymbol='%s'", private$target.gene)
            tbl <- dbGetQuery(db, query.string)
            new.order <- order(tbl$pvalue, decreasing=FALSE)
            tbl <- tbl[new.order,]
            dbDisconnect(db)
            private$tbl.eqtls <- tbl
            private$tbl.eqtls
            },

        getEQTLsInRegion = function(chrom, start, end){
            if(grepl("hagfish", Sys.info()[["nodename"]])){
               suppressWarnings(db.access.test <- try(system("/sbin/ping -c 1 khaleesi", intern=TRUE, ignore.stderr=TRUE)))
               if(length(db.access.test) == 0){
                  message(sprintf(("khaleesi unreachable, no eQTLS available")))
                  return(data.frame())
                  } # length 0
               } # on macos hagfish
            db <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="genereg2021", host="khaleesi")
            query.string <-  sprintf("select * from eqtls where chrom='%s' and hg38 > %d and hg38 < %d",
                               chrom, start, end)
            private$tbl.eqtls <- dbGetQuery(db, query.string)
            dbDisconnect(db)
            private$tbl.eqtls
            },

        getIdMap = function(){
           private$tbl.idMap
           },

        biomart.rsidToLoc = function(rsids){
           tbl.hg38 <- getBM(attributes=c("refsnp_id","allele","chr_name", "chrom_start"),
                             filters=c("snp_filter"),
                             values=list(rsids),
                             mart=private$snpMart.hg38)
           colnames(tbl.hg38) <- c("rsid", "allele", "chrom", "hg38")
           tbl.hg19 <- getBM(attributes=c("refsnp_id","allele","chr_name", "chrom_start"),
                             filters=c("snp_filter"),
                             values=list(rsids),
                             mart=private$snpMart.hg19)
           colnames(tbl.hg19) <- c("rsid", "allele", "chrom", "hg19")
           tbl.out <- merge(tbl.hg19[, c("chrom", "hg19", "rsid")],
                            tbl.hg38[, c("hg38", "rsid")], by="rsid")
           if(!grepl("chr", tbl.out$chrom[1]))
              tbl.out$chrom <- sprintf("chr%s", tbl.out$chrom)
           tbl.out[, c("chrom", "hg19", "hg38", "rsid")]
           }, # biomart.rsidToLoc

        rsidToLoc = function(rsids){
            # message(sprintf("--- EndophenotypeExplorer$rsidToLoc, starting time-consuming queries to SNPlocs"))
            rsids <- grep("^rs", rsids, value=TRUE)
            x.hg19 <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, rsids, ifnotfound="drop")
            tbl.hg19 <- as.data.frame(x.hg19) # do.call(rbind, lapply(x.hg19, as.data.frame))
            x.hg38 <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, rsids, ifnotfound="drop")
            tbl.hg38 <- as.data.frame(x.hg38) # do.call(rbind, lapply(x.hg38, as.data.frame))
            colnames(tbl.hg19)[1] <- "chrom"
            colnames(tbl.hg19)[2] <- "hg19"
            colnames(tbl.hg19)[4] <- "rsid"
            tbl.hg19$chrom <- as.character(tbl.hg19$chrom)
            colnames(tbl.hg38)[2] <- "hg38"
            colnames(tbl.hg38)[4] <- "rsid"
            tbl.out <- merge(tbl.hg19[, c("chrom", "hg19", "rsid")], tbl.hg38[, c("hg38", "rsid")], by="rsid")
            tbl.out[, c("chrom", "hg19", "hg38", "rsid")]
            }, # rsidToLoc

        new.rsidToLoc = function(rsids){
            # message(sprintf("--- EndophenotypeExplorer$rsidToLoc, starting time-consuming queries to SNPlocs"))
            rsids <- grep("^rs", rsids, value=TRUE)
            gr.hg19 <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, rsids, ifnotfound="drop")
            tbl.hg19 <- as.data.frame(gr.hg19)
            gr.hg38 <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh38, rsids, ifnotfound="drop")
            tbl.hg38 <- as.data.frame(gr.hg38)
            colnames(tbl.hg19)[1] <- "chrom"
            colnames(tbl.hg19)[2] <- "hg19"
            colnames(tbl.hg19)[4] <- "rsid"
            tbl.hg19$chrom <- as.character(tbl.hg19$chrom)
            colnames(tbl.hg38)[2] <- "hg38"
            colnames(tbl.hg38)[4] <- "rsid"
            tbl.out <- merge(tbl.hg19[, c("chrom", "hg19", "rsid")], tbl.hg38[, c("hg38", "rsid")], by="rsid")
            tbl.out[, c("chrom", "hg19", "hg38", "rsid")]
            }, # new.rsidToLoc

        biomart.locsToRSID = function(locs, genome){
              # expect locs in this form:  "2:127084193_G/A"
              # todo: enforce this, fail or accomodate if otherwise
            chroms <- unlist(lapply(strsplit(locs, ":"), "[", 1))
            loc.strings <- unlist(lapply(strsplit(locs, ":"), "[", 2))
            locs.base <- as.integer(sub("_.*$", "", loc.strings))
               # convert location string to the standard from required by biomart
               # which is chr:start:end
            locs.biomart <- sprintf("%s:%d:%d", chroms, locs.base, locs.base)
            if(genome == "hg19")
               tbl.out <-  getBM(attributes=c("refsnp_id","allele","chr_name", "chrom_start"),
                                 filters="chromosomal_region", # chr_name","start","end"),
                                 #values=locs.biomart,
                                 #values=c("2:127084193:127084193", "2:127084235:127084235"),
                                        #values=values,
                                 values=locs.biomart,
                                 mart=private$snpMart.hg19)
            if(genome == "hg38")
               tbl.out <-  getBM(attributes=c("refsnp_id","allele","chr_name", "chrom_start"),
                                 filters="chromosomal_region", # chr_name","start","end"),
                                 values=locs.biomart,
                                 #values=c("2:127084193:127084193", "2:127084235:127084235"),
                                 mart=private$snpMart.hg38)
            colnames(tbl.out) <- c("rsid", "allele", "chrom", genome)
            x <- tbl.out$rsid
            names(x) <- locs
            x
            }, # biomart.locsToRSID

        locsToRSID = function(locs, genome){
              # expect locs in this form:  "2:127084193_G/A"
              # todo: enforce this, fail or accomodate if otherwise

            regex <-  "(?<chrom>.*):(?<pos>\\d+)_(?<ref>.*)\\/(?<alt>.*)"
            matches <- as.integer(regexpr(regex, locs, perl=TRUE))
            good.locs <- which(matches == 1)
            bad.locs  <- which(matches == -1)
            result <- as.list(vector(mode="character", length=length(locs)))
            names(result) <- locs
            result[bad.locs] <- locs[bad.locs]
            good.loc.strings <- locs[good.locs]
            good.loc.strings <- sub("^chr", "", good.loc.strings)
            chroms <- unlist(lapply(strsplit(good.loc.strings, ":"), "[", 1))
            loc.strings <- unlist(lapply(strsplit(good.loc.strings, ":"), "[", 2))
            locs.base <- as.numeric(sub("_.*$", "", loc.strings))
            snplocs <- switch(genome,
                   "hg19" = SNPlocs.Hsapiens.dbSNP144.GRCh37,
                   "hg38" = SNPlocs.Hsapiens.dbSNP155.GRCh38
                   )
            gr <- GRanges(seqnames=chroms, IRanges(start=locs.base, end=locs.base))
            gr.snps <- snpsByOverlaps(snplocs, gr)
            tbl.rsids <- as.data.frame(gr.snps)[, c("seqnames", "pos", "RefSNP_id")]
            colnames(tbl.rsids)[1] <- "chrom"
            colnames(tbl.rsids)[2] <- "loc"
            colnames(tbl.rsids)[3] <- "rsid"
            tbl.rsids$chrom <- as.character(tbl.rsids$chrom)
            tbl.rsids$signature <- sprintf("%s:%s", tbl.rsids$chrom, tbl.rsids$loc)
            sig <- sprintf("%s:%s", chroms, locs.base)
            tbl.all <- data.frame(chrom=chroms, loc=locs.base,
                                  signature=sig,
                                  stringsAsFactors=FALSE)
            tbl.new <- merge(tbl.all, tbl.rsids[, c("rsid", "signature")], by="signature", all.x=TRUE)
            result[good.locs] <- tbl.new$rsid
            failed.lookups <- as.integer(which(is.na(result)))
            if(length(failed.lookups) > 0)
               result[failed.lookups] <- names(result)[failed.lookups]
            result
            }, # locsToRSID

        mapSampleIdToPatientAndCohort = function(sampleID){
           subset(private$tbl.idMap, sample==sampleID)
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
            tbl.c.mayo <- self$standardizeMayoPatientTable(private$tbl.clinical.mayo)
            tbl.c.rosmap <- self$standardizeRosmapPatientTable(private$tbl.clinical.rosmap)
            tbl.c.sinai <- self$standardizeSinaiPatientTable(private$tbl.clinical.sinai)
            stopifnot(all(colnames(tbl.c.mayo) == colnames(tbl.c.rosmap)))
            stopifnot(all(colnames(tbl.c.mayo) == colnames(tbl.c.sinai)))

            private$tbl.clinical <- rbind(tbl.c.mayo, tbl.c.rosmap, tbl.c.sinai)

            dir <- system.file(package="EndophenotypeExplorer", "extdata", "idMapping")
               # the crucial sample-to-patient mapping
               # this was difficult to obtain.  see the mapToPatientID function in
               # ~/github/TrenaProjectAD/explore/ampad.eQTLS/ldlr.R
               # todo: move this code to a prep directory in this package
            if(private$verbose) message(sprintf("--- about to load id mapping file"))
            if(private$verbose) message(sprintf("    dir: %s", dir))
            full.path <- file.path(dir, "tbl.sampleToPatientMap-rosmap-mayo-sinai.RData")
            #full.path <- file.path(dir, "tbl.vcfToPatientIDs.RData")
            if(private$verbose) message(sprintf("    full.path: %s", full.path))
            private$tbl.idMap <- get(load(full.path))
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

        sampleID.to.clinicalTable = function(sampleID){
            tbl.sub <- subset(private$tbl.idMap,  sample==sampleID)
            if(nrow(tbl.sub) == 0)
                return(data.frame())
            patientID <- tbl.sub$patient[1]
            study <- tbl.sub$study[1]
            tbl.patient <- switch(study,
              "mayo"   = self$standardizeMayoPatientTable(subset(private$tbl.clinical.mayo,
                                                          individualID==patientID)),
              "rosmap" = self$standardizeRosmapPatientTable(subset(private$tbl.clinical.rosmap,
                                                            individualID==patientID)),
              "sinai"  = self$standardizeSinaiPatientTable(subset(private$tbl.clinical.sinai,
                                                           individualID==patientID))
              )
            if(nrow(tbl.patient) == 0) return(data.frame())
            tbl.patient$study <- study
            tbl.patient$sampleID <- sampleID
            coi <- c("patientID", "sampleID", "study", "study","sex","ethnicity","apoeGenotype",
                     "braak","cerad","cogdx","pmi", "ageAtDeath")
            tbl.patient[, coi]
            },

        getStandardClinicalColumnNames = function(){
            private$standard.clinical.columns
            },

        standardizeMayoPatientTable = function(tbl){
           tbl$ageDeath <- as.numeric(sub("+", "", tbl$ageDeath, fixed=TRUE))

           coi <- c("individualID", "individualIdSource", "sex", "ethnicity", "apoeGenotype",
                    "Braak", "CERAD", "pmi", "ageDeath")
           tbl.0 <- tbl[, coi]
           tbl.0$cogdx <- NA
           colnames(tbl.0) <- private$standard.clinical.columns
           tbl.0
           },

        standardizeRosmapPatientTable = function(tbl){
           tbl$age_death <- as.numeric(sub("+", "", tbl$age_death, fixed=TRUE))
           coi <- c("individualID", "Study", "msex", "race", "apoe_genotype", "braaksc", "ceradsc",
                    "pmi", "age_death", "cogdx")
           private$standard.clinical.columns <- c("patientID", "study", "sex","ethnicity",
                                                  "apoeGenotype","braak","cerad", "pmi",
                                                  "ageAtDeath", "cogdx")


           tbl.0 <- tbl[, coi]
           colnames(tbl.0) <- private$standard.clinical.columns
           tbl.0
           },

        standardizeSinaiPatientTable = function(tbl){
           tbl$ageDeath <- as.numeric(sub("+", "", tbl$ageDeath, fixed=TRUE))
           coi <- c("individualID","individualIdSource","sex","ethnicity","apoeGenotype",
                    "Braak","CERAD","pmi","ageDeath")
           tbl.0 <- tbl[, coi]
           tbl.0$cogdx <- NA
           colnames(tbl.0) <- private$standard.clinical.columns
           tbl.0
           },

        getAggregatedAlleleFrequencies = function(rsid, quiet=TRUE){

            uri <- sprintf("https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/%s/frequency",
                  sub("rs", "", rsid))
            x <- NULL
            tries <- 0
            max.tries <- 5
            while(is.null(x) & tries <= max.tries){
              tries <- tries + 1
              if(!quiet)
                 message(sprintf("requesting allele frequencies for %s, iteration %d", rsid, tries))
              response <- GET(uri)
              suppressMessages(x <- fromJSON(httr::content(response, as="text"))$results)
              }
            if(is.null(x)){
                message(sprintf("failed to get allele frequencies for %s, %d tries",
                                rsid, tries))
                return(data.frame())
                }
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
           },

             # mayoMatrixTissueSuffix is either _TCX or _CER
             # use it to get only one sub matrix or the other, then remove it
        splitExpressionMatrixByMutationStatusAtRSID = function(mtx, rsid, study.name){
           stopifnot(study.name %in% c("mayo", "sinai", "rosmap"))
           result <- switch(study.name,
                             "mayo" = self$splitMayoRnaMatrixByGenotype(mtx,rsid),
                             "sinai" = self$splitSinaiRnaMatrixByGenotype(mtx,rsid),
                             "rosmap" = self$splitRosmapRnaMatrixByGenotype(mtx,rsid))
           return(result)
           }, # splitExpressionMatrixByMutationStatusAtRSID

        #------------------------------------------------------------
        splitMayoRnaMatrixByGenotype = function(mtx,rsid){
           tbl.map <- self$getIdMap()
           mtx.geno <- self$getGenoMatrixByRSID(rsid)
           dim(mtx.geno)
           table(mtx.geno)     # 794  861 239   for rs28834970
           samples.hom <- names(which(mtx.geno[1,] == "1/1"))
           samples.het <- names(which(mtx.geno[1,] == "0/1"))
           samples.wt  <- names(which(mtx.geno[1,] == "0/0"))

           patients.wt  <- unique(subset(tbl.map, sample %in% samples.wt &
                                                  study=="mayo" &
                                                  assay=="vcf")$patient)
           patients.hom <- unique(subset(tbl.map, sample %in% samples.hom &
                                                  study=="mayo" &
                                                  assay=="vcf")$patient)
           patients.het <- unique(subset(tbl.map, sample %in% samples.het &
                                                  study=="mayo" &
                                                  assay=="vcf")$patient)
           patients.mut <- unique(sort(unique(c(patients.hom, patients.het))))

           rna.samples.wt <- subset(tbl.map, patient %in% patients.wt &
                                             study=="mayo" &
                                             assay=="rnaseq")$sample
           rna.samples.hom <- subset(tbl.map, patient %in% patients.hom &
                                              study=="mayo" &
                                              assay=="rnaseq")$sample
           rna.samples.het <- subset(tbl.map, patient %in% patients.het &
                                              study=="mayo" &
                                              assay=="rnaseq")$sample
           rna.samples.mut <- subset(tbl.map, patient %in% patients.mut &
                                              study=="mayo" &
                                              assay=="rnaseq")$sample
           rna.samples.wt <- sub("_TCX", "", rna.samples.wt)
           rna.samples.wt <- sub("_CER", "", rna.samples.wt)
           rna.samples.mut <- sub("_TCX", "", rna.samples.mut)
           rna.samples.mut <- sub("_CER", "", rna.samples.mut)
           rna.samples.het <- sub("_TCX", "", rna.samples.het)
           rna.samples.het <- sub("_CER", "", rna.samples.het)
           rna.samples.hom <- sub("_TCX", "", rna.samples.hom)
           rna.samples.hom <- sub("_CER", "", rna.samples.hom)
           rna.samples.wt  <- intersect(colnames(mtx), rna.samples.wt)
           rna.samples.mut <- intersect(colnames(mtx), rna.samples.mut)
           rna.samples.hom <- intersect(colnames(mtx), rna.samples.hom)
           rna.samples.het <- intersect(colnames(mtx), rna.samples.het)
           mtx.hom <- mtx[, rna.samples.hom]
           mtx.het <- mtx[, rna.samples.het]
           mtx.wt <-  mtx[, rna.samples.wt]
           mtx.mut <- mtx[, rna.samples.mut]
           patient.distribution <- list(wt=ncol(mtx.wt),
                                        mut=ncol(mtx.mut),
                                        het=ncol(mtx.het),
                                        hom=ncol(mtx.hom))
           return(list(wt=mtx.wt, mut=mtx.mut, het=mtx.het, hom=mtx.hom,
                       genotypes.rna=patient.distribution,
                       genotypes.vcf=list(wt=length(patients.wt),
                                          mut=length(patients.mut),
                                          het=length(patients.het),
                                          hom=length(patients.hom))
                       ))
           }, # splitMayo

        #------------------------------------------------------------
        splitSinaiRnaMatrixByGenotype = function(mtx,rsid){
           tbl.map <- self$getIdMap()
           mtx.geno <- self$getGenoMatrixByRSID(rsid)
           dim(mtx.geno)
           table(mtx.geno)     # 794  861 239   for rs28834970
           samples.hom <- names(which(mtx.geno[1,] == "1/1"))
           samples.het <- names(which(mtx.geno[1,] == "0/1"))
           samples.wt  <- names(which(mtx.geno[1,] == "0/0"))

           patients.wt  <- unique(subset(tbl.map,
                                         sample %in% samples.wt &
                                         study=="sinai" &
                                         assay=="vcf")$patient)
           patients.hom  <- unique(subset(tbl.map,
                                         sample %in% samples.hom &
                                         study=="sinai" &
                                         assay=="vcf")$patient)
           patients.het <- unique(subset(tbl.map,
                                         sample %in% samples.het &
                                         study=="sinai" &
                                         assay=="vcf")$patient)
           patients.mut <- unique(sort(unique(c(patients.hom, patients.het))))

           rna.samples.wt  <- intersect(colnames(mtx), patients.wt)
           rna.samples.mut <- intersect(colnames(mtx), patients.mut)
           rna.samples.hom <- intersect(colnames(mtx), patients.hom)
           rna.samples.het <- intersect(colnames(mtx), patients.het)
           mtx.hom <- mtx[, rna.samples.hom]
           mtx.het <- mtx[, rna.samples.het]
           mtx.wt <-  mtx[, rna.samples.wt]
           mtx.mut <- mtx[, rna.samples.mut]
           patient.distribution <- list(wt=ncol(mtx.wt),
                                        mut=ncol(mtx.mut),
                                        het=ncol(mtx.het),
                                        hom=ncol(mtx.hom))

           return(list(wt=mtx.wt, mut=mtx.mut, het=mtx.het, hom=mtx.hom,
                       genotypes.rna=patient.distribution,
                       genotypes.vcf=list(wt=length(patients.wt),
                                          mut=length(patients.mut),
                                          het=length(patients.het),
                                          hom=length(patients.hom))
                       ))
           }, # splitSinaiRnaMatrixByGenotype

        #------------------------------------------------------------
        splitRosmapRnaMatrixByGenotype = function(mtx, rsid){
           tbl.map <- self$getIdMap()
           mtx.geno <- self$getGenoMatrixByRSID(rsid)
           dim(mtx.geno)
           table(mtx.geno)     # 794  861 239   for rs28834970
           samples.hom <- names(which(mtx.geno[1,] == "1/1"))
           samples.het <- names(which(mtx.geno[1,] == "0/1"))
           samples.wt  <- names(which(mtx.geno[1,] == "0/0"))

           patients.wt  <- unique(subset(tbl.map,
                                         sample %in% samples.wt &
                                         study=="rosmap" &
                                         assay=="vcf")$patient)
           patients.hom  <- unique(subset(tbl.map,
                                         sample %in% samples.hom &
                                         study=="rosmap" &
                                         assay=="vcf")$patient)
           patients.het <- unique(subset(tbl.map,
                                         sample %in% samples.het &
                                         study=="rosmap" &
                                         assay=="vcf")$patient)
           patients.mut <- unique(sort(unique(c(patients.hom, patients.het))))

           rna.samples.wt  <- intersect(colnames(mtx), patients.wt)
           rna.samples.mut <- intersect(colnames(mtx), patients.mut)
           rna.samples.hom <- intersect(colnames(mtx), patients.hom)
           rna.samples.het <- intersect(colnames(mtx), patients.het)
           mtx.hom <- mtx[, rna.samples.hom]
           mtx.het <- mtx[, rna.samples.het]
           mtx.wt <-  mtx[, rna.samples.wt]
           mtx.mut <- mtx[, rna.samples.mut]
           patient.distribution <- list(wt=ncol(mtx.wt),
                                        mut=ncol(mtx.mut),
                                        het=ncol(mtx.het),
                                        hom=ncol(mtx.hom))

           return(list(wt=mtx.wt, mut=mtx.mut, het=mtx.het, hom=mtx.hom,
                       genotypes.rna=patient.distribution,
                       genotypes.vcf=list(wt=length(patients.wt),
                                          mut=length(patients.mut),
                                          het=length(patients.het),
                                          hom=length(patients.hom))
                       ))

            }, # splitRosmapRnaMatrixByGenotype

        #------------------------------------------------------------
        trenaScoreGenotypeStratifiedExpression = function(mtx.1, mtx.2, targetGene, tfs){
           mtx.1[is.na(mtx.1)] <- 0
           mtx.2[is.na(mtx.2)] <- 0
           solver.names <- c("spearman", "pearson", "bicor", "randomForest", "xgboost",
                             "lasso", "ridge")
           solver <- EnsembleSolver(mtx.1,
                                    targetGene=targetGene,
                                    candidateRegulators=tfs,
                                    solverNames=solver.names,
                                    geneCutoff=0.9)
           tbl.trena.1 <- run(solver)
           tbl.trena.1 <- tbl.trena.1[order(abs(tbl.trena.1$spearman), decreasing=TRUE),]

           solver <- EnsembleSolver(mtx.2,
                                    targetGene=targetGene,
                                    candidateRegulators=tfs,
                                    solverNames=solver.names,
                                    geneCutoff=0.9)
           tbl.trena.2 <- run(solver)
           tbl.trena.2 <- tbl.trena.2[order(abs(tbl.trena.2$spearman), decreasing=TRUE),]

           calculateModelDeltas <- function(tbl.trena.1, tbl.trena.2, scoreName){
               tbl.delta <- data.frame(tf=sort(unique(c(tbl.trena.1$gene, tbl.trena.2$gene))),
                                       method=scoreName, stringsAsFactors=FALSE)
               tbl.delta$wt <- 0
               tbl.delta$mut <- 0

               indices <- match(tbl.trena.1$gene, tbl.delta$tf)
               tbl.delta$wt[indices] <- tbl.trena.1[, scoreName]
               indices <- match(tbl.trena.2$gene, tbl.delta$tf)
               tbl.delta$mut[indices] <- tbl.trena.2[, scoreName]
               tbl.delta[, "delta"] <- tbl.delta$mut - tbl.delta$wt
               new.order <- order(abs(tbl.delta$delta), decreasing=TRUE)
               tbl.delta <- tbl.delta[new.order,]
               rownames(tbl.delta) <- NULL
               tbl.delta
               } # calculateModelDeltas

           tbl.bicor <- calculateModelDeltas(tbl.trena.1, tbl.trena.2, "bicor")
           tbl.rf    <- calculateModelDeltas(tbl.trena.1, tbl.trena.2, "rfScore")
           tbl.spear <- calculateModelDeltas(tbl.trena.1, tbl.trena.2, "spearmanCoeff")
           tbl.lasso <- calculateModelDeltas(tbl.trena.1, tbl.trena.2, "betaLasso")
           tbl.ridge <- calculateModelDeltas(tbl.trena.1, tbl.trena.2, "betaRidge")
           tbl.pearson <- calculateModelDeltas(tbl.trena.1, tbl.trena.2, "pearsonCoeff")
           tbl.xgboost <- calculateModelDeltas(tbl.trena.1, tbl.trena.2, "xgboost")
           list(trena.1=tbl.trena.1, trena.2=tbl.trena.2,
                bicor=tbl.bicor, rfScore=tbl.rf,
                spearmanCoeff=tbl.spear, pearsonCoeff=tbl.pearson,
                betaLasso=tbl.lasso, betaRidge=tbl.ridge, xgboost=tbl.xgboost)
           }, # trenaScoreGenotypeStratifiedExpression


           # duplicates parts of splitExpressionMatrixByMutationStatusAtRSID
        subsetAndRelabelGenoMatrixByPatientIDs = function(mtx.geno, studyName){
           tbl.map <- self$getIdMap()
           stopifnot(studyName %in% tbl.map$study)
           tbl.sub <- subset(tbl.map,
                             study==studyName &
                             sample %in% colnames(mtx.geno) &
                             assay=="vcf")
           if(nrow(tbl.sub) == 0)
               return(NA)
           mtx.geno.sub <- mtx.geno[, tbl.sub$sample, drop=FALSE]
           indices <- match(colnames(mtx.geno.sub), tbl.sub$sample)
           colnames(mtx.geno.sub) <- tbl.sub$patient[indices]
           mtx.geno.sub
           }, # subsetAndRelabelGenoMatrixByPatientIDs

           # fisher's exact and T test on enrichment of braak >= 5, <=1
           # in rosmpa patients, when separated on having/not-having the
           # specified variant
        extremeBraakScoreSignificanceAtRSID = function(rsid){
            result <- list(pval.t=NA,
                           pval.fisher=NA,
                           tbl.geno=data.frame(),
                           tbl.pt=data.frame(),
                           mtx.geno=matrix(),
                           mtx.geno.study=matrix(),
                           pt.ad=c(),
                           pt.ctl=c())

            mtx.geno.1 <- self$getGenoMatrixByRSID(rsid)
            if(all(is.na(mtx.geno.1)))
                return(result)

            new.names <- etx$locsToRSID(rownames(mtx.geno.1), "hg19")
            rownames(mtx.geno.1) <- new.names
            mtx.geno.pt.rosmap <- self$subsetAndRelabelGenoMatrixByPatientIDs(mtx.geno.1, "rosmap")
            if(all(is.na(mtx.geno.pt.rosmap)))
                return(result)

            tbl.pt <- self$get.rosmap.patient.table(NA)
            rosmap.patients <- intersect(tbl.pt$individualID, colnames(mtx.geno.pt.rosmap))
            length(rosmap.patients)   # 1143

            tbl.pt.rosmap <- subset(tbl.pt, individualID %in% rosmap.patients)
            dim(tbl.pt.rosmap)  # 1143 18

            dim(mtx.geno.pt.rosmap)
            table(mtx.geno.pt.rosmap)  # mtx.geno.pt.rosmap
                                        # 0/0 0/1 1/1
                                        # 827 296  28
            pt.ad <-  subset(tbl.pt.rosmap, braaksc >=5)$individualID   # 303
            pt.ctl <- subset(tbl.pt.rosmap, braaksc <=1)$individualID  # 85


            table(mtx.geno.pt.rosmap)  # mtx.geno.pt.rosmap
                                        #  0/0 0/1 1/1
                                        #  827 296  28

            pt.ad <-  subset(tbl.pt.rosmap, braaksc >=5)$individualID   # 303
            pt.ctl <- subset(tbl.pt.rosmap, braaksc <=1)$individualID  # 85

            ad.wt <- 0
            at.het <- 0
            ad.hom <- 0

            ad.wt  <- length(grep("0/0", (mtx.geno.pt.rosmap[, pt.ad])))
            ad.het <- length(grep("0/1", (mtx.geno.pt.rosmap[, pt.ad])))
            ad.hom <- length(grep("1/1", (mtx.geno.pt.rosmap[, pt.ad])))

            ctl.wt  <- length(grep("0/0", (mtx.geno.pt.rosmap[, pt.ctl])))
            ctl.het <- length(grep("0/1", (mtx.geno.pt.rosmap[, pt.ctl])))
            ctl.hom <- length(grep("1/1", (mtx.geno.pt.rosmap[, pt.ctl])))

            tbl.summary <- data.frame(wt=c(ad.wt, ctl.wt), het=c(ad.het, ctl.het), hom=c(ad.hom, ctl.hom),
                                      row.names=c("ad", "ctl"))
            print(tbl.summary)

            ad.vector <- with(tbl.summary["ad",],
                              c(rep(0, wt),
                                rep(1, het),
                                rep(2, hom)))
            ctl.vector <- with(tbl.summary["ctl",],
                               c(rep(0, wt),
                                 rep(1, het),
                                 rep(2, hom)))

            pval.t <- tryCatch({
                t.test(ad.vector, ctl.vector)$p.value
                }, error=function(e){return(1)})

            pval.fisher <- tryCatch({
                fisher.test(tbl.summary)$p.value
                }, error=function(e){return(1)})

            return(list(pval.t=pval.t,
                        pval.fisher=pval.fisher,
                        tbl.geno=tbl.summary,
                        tbl.pt=tbl.pt.rosmap,
                        mtx.geno=mtx.geno.1,
                        mtx.geno.study=mtx.geno.pt.rosmap,
                        pt.ad=pt.ad,
                        pt.ctl=pt.ctl))

            } # extremeBraakScoreSignificanceAtRSID

        #------------------------------------------------------------
        ) # public

    ) # class EndophenotypeExplorer

#--------------------------------------------------------------------------------


