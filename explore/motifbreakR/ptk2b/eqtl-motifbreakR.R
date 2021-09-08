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
rsids.all <- tbl.eQTL$rsid
length(rsids.all)
head(tbl.eQTL)
dim(tbl.eQTL)
dim(subset(tbl.eQTL, pvalue < 0.005)) 
tbl.eQTL[1:5,]
subset(tbl.eQTL, pvalue < 0.005)

direction <- "bottom"
#direction <- "random"
#direction <- "middle"
#direction <- "top"

#rich.pvals.end <- nrow(subset(tbl.eQTL, pvalue <= 0.05))
#rsids.oi <- rsids.all[1000:2000]
#pval.threshold <- 0.005
#rsids.oi <- subset(tbl.eQTL, pvalue <= pval.threshold)$rsid
rsids.oi <- tail(tbl.eQTL$rsid, n=200)
match(rsids.oi, tbl.eQTL$rsid)

count <- length(rsids.oi)
printf("found %d rsids", count)



#rsids.midsig <- subset(tbl.eQTL, pvalue <= 0.55 & pvalue >= 0.50)$rsid
#length(rsids.midsig)
#rsids.midsig <- sample(rsids.midsig, size=40)

#snps.gr <- snps.from.rsid(rsid=rsids.oi,
#snps.gr <- snps.from.rsid(rsid=rsids.sig,
printf("--- getting snps")
print(head(rsids.oi))
snps.gr <- snps.from.rsid(rsid=rsids.oi, 
                          dbSNP=SNPlocs.Hsapiens.dbSNP151.GRCh38,
                          search.genome=BSgenome.Hsapiens.UCSC.hg38)

printf("----- got snps: %d", length(snps.gr))
length(snps.gr)
printf("----- get motifs")

#motifs.selected <- get(load("motifs.50.80.pkt2b.RData"))
#motifs.selected <- get(load("motifs.38tfs.60motifs.pkt2b.RData"))
motifs.selected <- get(load("motifs.2tfs.3motifs.pkt2b.RData"))

length(motifs.selected)

worker.count <- min(length(rsids.oi), 40)
printf("worker.count: %d", worker.count)
#param <- SnowParam(workers=worker.count, type = "SOCK")
bpparam <- MulticoreParam(workers=worker.count)


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
print(colSums(tbl.summary))

filename <- sprintf("ptk2b.eqtls.%s.%d.%s.motifbreaks.RData", direction, count,
                    format(Sys.time(), "%d%b.%H%M"))
save(results, tbl.breaks, tbl.summary, file=filename)

# 
#  rsids.to.use <- rsids.all[sort(sample(seq_len(length(rsids)), size=3))]
#  
#  run.breaker <- function(snp.gr){
#      printf("run.breaker on %s", paste(rsids, collapse= " "))
#      results <- motifbreakR(snpList = snp.gr,
#                             filterp = TRUE,
#                             pwmList = motifs.selected,
#                             show.neutral=FALSE,
#                             method = c("ic", "log", "notrans")[1],
#                             bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
#                             verbose=TRUE)
#      }
#  
#  direction <- "random"
#  count <- 5
#  if(direction == "top")
#    rsids <- head(rsids.all, n=count)
#  if(direction == "bottom")
#    rsids <- tail(rsids.all, n=count)
#  if(direction == "random")
#    rsids <- rsids.all[sort(sample(seq_len(length(rsids.all)), size=count))]
#  
#  snps.gr <- snps.from.rsid(rsid=rsids,
#                            dbSNP=SNPlocs.Hsapiens.dbSNP151.GRCh38,
#                            search.genome=BSgenome.Hsapiens.UCSC.hg38)
#  x <- bplapply(list(snps.gr[[1]], snps.gr[[2]]), run.breaker, BPParame=MulitcoreParam()))
#  #x <- bplapply(rsids, run.breaker, BPParame=MulitcoreParam()))
#  
