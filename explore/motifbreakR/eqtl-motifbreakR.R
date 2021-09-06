library(EndophenotypeExplorer)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifbreakR)
library(MotifDb)

etx <- EndophenotypeExplorer$new("PTK2B", "hg38");
tbl.eQTL <- etx$getEQTLsForGene()

rsids <- subset(tbl.eQTL, pvalue < 0.001)$rsid
length(rsids)  # 17


snps.gr <- snps.from.rsid(rsid=rsids,
                          dbSNP=SNPlocs.Hsapiens.dbSNP151.GRCh38,
                          search.genome=BSgenome.Hsapiens.UCSC.hg38)


motifs.selected <- query(MotifDb, "sapiens", c("jaspar2018", "hocomoco-core-A"))

results <- motifbreakR(snpList = snps.gr,
                       filterp = TRUE,
                       pwmList = motifs.selected,
                       show.neutral=FALSE,
                       method = c("ic", "log", "notrans")[1],
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam(),
                       verbose=TRUE)

tbl <- as.data.frame(results, row.names=NULL)
tbl$pctDelta <- tbl$pctAlt - tbl$pctRef
tbl$scoreDelta <- tbl$scoreAlt - tbl$scoreRef
cor(tbl$pctDelta, tbl$scoreDelta)   # 0.94
fivenum(abs(tbl$scoreDelta))        # 0.4001632 0.6700027 1.0379065 1.5177660 1.9998288
tbl.strongest <- subset(tbl, scoreDelta > 1.995)
dim(tbl.strongest)



