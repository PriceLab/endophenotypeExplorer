library(EndophenotypeExplorer)
library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifbreakR)
library(MotifDb)
library(BiocParallel)

etx <- EndophenotypeExplorer$new("PTK2B", "hg38");
tbl.eQTL <- etx$getEQTLsForGene()
tbl.eQTL  <- subset(tbl.eQTL, study=="ampad-mayo" & tissue=="tcx")
new.order <- order(tbl.eQTL$pvalue, decreasing=FALSE)
tbl.eQTL <- tbl.eQTL[new.order,]
head(tbl.eQTL)
dim(tbl.eQTL)

#selected <- sort(sample(seq_len(nrow(tbl.eQTL)), size=100))
# tbl.eQTL.small <- tbl.eQTL[selected,]
# dim(tbl.eQTL.small)
# rsids <- unique(tbl.eQTL.small$rsid)
rsids <- head(tbl.eQTL$rsid, n=10)
subset(tbl.eQTL, rsid %in% rsids)

length(rsids)  # 17


printf("----- get snps")
snps.gr <- snps.from.rsid(rsid=rsids,
                          dbSNP=SNPlocs.Hsapiens.dbSNP151.GRCh38,
                          search.genome=BSgenome.Hsapiens.UCSC.hg38)
length(snps.gr)
printf("----- get motifs")

motifs.selected <- get(load("motifs.50.80.pkt2b.RData"))

worker.count <- min(length(rsids), 40)
#param <- SnowParam(workers=worker.count, type = "SOCK")
bpparam <- MulticoreParam()


printf("----- running breakR")
results <- motifbreakR(snpList = snps.gr,
                       filterp = TRUE,
                       pwmList = motifs.selected,
                       show.neutral=FALSE,
                       method = c("ic", "log", "notrans")[1],
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = bpparam,
                       verbose=TRUE)
tbl.breaks <- as.data.frame(results, row.names=NULL)
dim(tbl.breaks)
tbl.breaks <- subset(tbl.breaks, effect=="strong")
dim(tbl.breaks)
tbl.breaks$pctDelta <- tbl.breaks$pctAlt - tbl.breaks$pctRef
tbl.breaks$scoreDelta <- tbl.breaks$scoreAlt - tbl.breaks$scoreRef
xtab.summary <- with(tbl.breaks, table(geneSymbol, SNP_id))
tbl.summary <- as.data.frame.matrix(xtab.summary)
colSums(tbl.summary)

                         






rownames(tbl.summary) <- tbl.summary$geneSymbol
printf("----- breakR complete")

save(results, file=sprintf("ptk2b.eqtls.%d.motifbreaks.RData", length(rsids)))
