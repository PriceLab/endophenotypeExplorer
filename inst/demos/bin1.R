library(VariantAnnotation)
library(EndophenotypeExplorer)
library(RUnit)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(MafDb.1Kgenomes.phase3.GRCh38)
library(MafDb.1Kgenomes.phase3.hs37d5)

vcf.base.url <- "https://igv-data.systemsbiology.net/static/ampad"
vcf.data.file.chr2 <- "NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_2.recalibrated_variants.vcf.gz"
vcf.url <- sprintf("%s/%s", vcf.base.url, vcf.data.file.chr2)


etx <- EndophenotypeExplorer$new("BIN1", "hg19", vcf.url)

# ATAC-seq region at the 5' UTR of at least one (perhaps the primary) transcript.
# look at it first in hg38, using data from Xue at Mayo Jacksonville
tbl.atac <- get(load("~/github/TrenaProjectAD/explore/mayo-epigenetics/atac/dbaConsensusRegionsScored.74273x30.RData"))

igv <- start.igv("BIN1", "hg38")
roi <- getGenomicRegion(igv)

tbl.bin1 <- subset(tbl.atac, chrom==roi$chrom & start >= roi$start & end <= roi$end)
dim(tbl.bin1) # 8 31
track <- DataFrameAnnotationTrack("atac", tbl.bin1[, c("chrom", "start", "end")], color="red",
                                  trackHeight=25)
displayTrack(igv, track)
  # zoom into the 824 bp region which diffbind shows is differentially open in PSP and AD vs control
  # "chr2:127,084,077-127,084,907"

roi <- getGenomicRegion(igv)

eqtls.file <- file.path("~/github/TrenaProjectAD/explore/mayo-epigenetics/progressReport-june-2021",
                        "tbl.eqtls.all38genes.RData")
file.exists(eqtls.file)
tbl.eqtls.all <- get(load(eqtls.file))
tbl.eqtls.sub <- subset(tbl.eqtls.all, chrom==roi$chrom & hg38 >= roi$start & hg38 <= roi$end)
dim(tbl.eqtls.sub) # 15 11
tbl.eqtls.hg19 <- tbl.eqtls.sub[, c("chrom", "hg19", "hg19", "pvalue")]

tbl.track <- tbl.eqtls.sub[, c("chrom", "hg38", "hg38", "pvalue")]
tbl.track$pvalue <- -log10(tbl.track$pvalue)
tbl.track$hg38 <- tbl.track$hg38 - 1


track <- DataFrameQuantitativeTrack("eqtls-scored", tbl.track, color="blue", autoscale=TRUE)
displayTrack(igv, track)

tbl.track <- tbl.eqtls.sub[, c("chrom", "hg38", "hg38", "rsid")]
tbl.track$hg38 <- tbl.track$hg38 - 1
track <- DataFrameAnnotationTrack("eqtl-rsids", tbl.track, color="blue", trackHeigh=25)
displayTrack(igv, track)

subset(tbl.eqtls.sub, rsid %in% c("rs10200967", "rs17014923"))
 #        chrom      hg19      hg38       rsid       pvalue            ensg genesymbol        study tissue   assay     score
 # 8433    chr2 127841930 127084354 rs17014923 0.4269389599 ENSG00000136717       BIN1   ampad-mayo    tcx unknown 0.3696342
 # 116215  chr2 127841769 127084193 rs10200967 0.0008776247 ENSG00000136717       BIN1   ampad-mayo    cer unknown 3.0566912
 # 155819  chr2 127841930 127084354 rs17014923 0.0974907406 ENSG00000136717       BIN1   ampad-mayo    cer unknown 1.0110366
 # 263496  chr2 127841769 127084193 rs10200967 0.0062688610 ENSG00000136717       BIN1 ampad-rosmap  dlpfc unknown 2.2028114
 # 296407  chr2 127841930 127084354 rs17014923 0.0008052473 ENSG00000136717       BIN1 ampad-rosmap  dlpfc unknown 3.0940707
 # 393521  chr2 127841769 127084193 rs10200967 0.2376787731 ENSG00000136717       BIN1   ampad-mayo    tcx unknown 0.6240096

 # chrom      hg19      hg38       rsid       pvalue            ensg genesymbol        study tissue   assay     score
 #  chr2 127841930 127084354 rs17014923 0.4269389599 ENSG00000136717       BIN1   ampad-mayo    tcx unknown 0.3696342
 #  chr2 127841930 127084354 rs17014923 0.0974907406 ENSG00000136717       BIN1   ampad-mayo    cer unknown 1.0110366
 #  chr2 127841930 127084354 rs17014923 0.0008052473 ENSG00000136717       BIN1 ampad-rosmap  dlpfc unknown 3.0940707




library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37, "8:27220295-27220330")



roi <- GRanges(seqnames="2", IRanges(start=127841769, end=127841930))
x <- readVcf(vcf.url, "hg19", roi)
dim(geno(x)$GT)
mtx.geno <- geno(x)$GT
mtx.geno[, 1:10]
variants <- rownames(mtx.geno)
x <- lapply(variants, function(snp) data.frame(variant=snp, table(mtx.geno[snp,])))
do.call(rbind, x)

gr <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, c("rs17014923", "rs10200967"))

mafdb <- MafDb.1Kgenomes.phase3.hs37d5
gscores(mafdb, gr)
#   [1]        2 127841930      * |  rs17014923                Y      0.21
#   [2]        2 127841769      * |  rs10200967                Y      0.30

IUPAC_CODE_MAP[["Y"]]  # CT
