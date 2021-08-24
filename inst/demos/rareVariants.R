library(EndophenotypeExplorer)
library(plyr)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(ghdb)
ghdb <- GeneHancerDB()

tbl.williams <- get(load(system.file(package="EndophenotypeExplorer", "extdata", "gwas",
                                     "williams-natureNeuroscience2020.RData")))
dim(tbl.williams)
rsids <- unique(grep("^rs", unlist(strsplit(tbl.williams$rsid, ",")), value=TRUE))

tbl.posthuma <- get(load(system.file(package="EndophenotypeExplorer", "extdata", "gwas",
                                     "tbl.posthuma-38-geneAssociations-curated-3828x12.RData")))
dim(tbl.posthuma)

etx <- EndophenotypeExplorer$new("CNTNAP2", "hg19")
tbls <- list()
   for(rsid in rsids){
      message(sprintf("--- %s", rsid))
      tbl <- etx$getAggregatedAlleleFrequencies(rsid)
      tbls[[rsid]] <- tbl
      Sys.sleep(1)
      }
tbl <- do.call(rbind.fill, tbls)
dim(tbl)
tbl.rareEuropean <- subset(tbl, population=="European" & min.freq < 1)
tbl.rareEuropean
tbl.rareEuropean <- merge(tbl.rareEuropean, tbl.williams, by="rsid", all.x=TRUE)
new.order <- order(tbl.rareEuropean$min.freq, decreasing=FALSE)
tbl.rareEuropean <- tbl.rareEuropean[new.order,]

tbl.snpLocs <- as.data.frame(snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, tbl.rareEuropean$rsid))
colnames(tbl.snpLocs)[1] <- "chrom"
colnames(tbl.snpLocs)[4] <- "rsid"
colnames(tbl.snpLocs)[2] <- "hg19"
tbl.snpLocs$chrom <- as.character(tbl.snpLocs$chrom)

tbl.rareEuropean <- merge(tbl.rareEuropean, tbl.snpLocs[, c("chrom", "hg19", "rsid")], by="rsid")

snp.loc <- subset(tbl.rareEuropean, locusOrGene=="CNTNAP2")$hg19
shoulder <- 1000
mtx.geno <- etx$getGenoMatrix("7", snp.loc-shoulder, snp.loc+shoulder)
locs <- rownames(mtx.geno)
chroms <- unlist(lapply(strsplit(locs, ":"), "[", 1))
loc.strings <- unlist(lapply(strsplit(locs, ":"), "[", 2))
locs <- as.integer(sub("_.*$", "", loc.strings))

rsids.region <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, GRanges(seqnames=chroms, IRanges(start=locs, end=locs)))$RefSNP_id
length(rsids.region)

igv <- start.igv("CNTNAP2", "hg19")
tbl.all <- retrieveEnhancersFromDatabase(ghdb, "CNTNAP2", tissues="all")
dim(tbl.all) # 47 16
tbl.hg19 <- to.hg19(ghdb, tbl.all)
dim(tbl.hg19)

tbl.track <- tbl.hg19[, c("chrom", "start", "end", "combinedscore")]
tbl.track$combinedscore <- asinh(tbl.track$combinedscore)
track <- DataFrameQuantitativeTrack("gh", tbl.track, autoscale=TRUE, color="brown")
displayTrack(igv, track)

tbls.region <- list()
for(rsid in rsids.region){
    message(sprintf("--- %s", rsid))
    tbl <- etx$getAggregatedAlleleFrequencies(rsid)
    tbls.region[[rsid]] <- tbl
    Sys.sleep(1)
    }

tbl.region <- do.call(rbind.fill, tbls.region)
length(rsids)
dim(tbl.region)
dim(subset(tbl.region, population=="European" & min.freq < 1))

unlist(lapply(strsplit(locs, ":")

lapply(seq_len(nrow(mtx.geno)), function(r) table(mtx.geno[r,]))

dim(mtx.geno)
