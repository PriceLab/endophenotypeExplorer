library(EndophenotypeExplorer)

etx <- EndophenotypeExplorer$new("PTK2B", "hg38");
tbl.eQTL <- etx$getEQTLsForGene()
tbl.eQTL  <- subset(tbl.eQTL, study=="ampad-mayo" & tissue=="tcx")
new.order <- order(tbl.eQTL$pvalue, decreasing=FALSE)
tbl.eQTL <- tbl.eQTL[new.order,]

print(load("ptk2b.eqtls.20.motifbreaks.RData"))
print(load("ptk2b.eqtls.top.14.08Sep.0947.motifbreaks.RData"))
results.t20 <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.t20 <- tbl.breaks
tbl.summary.t20 <- tbl.summary
tbl.breaks.strong.t20  <- subset(tbl.breaks.t20, pctRef > 0.9)
sort(table(tbl.breaks.strong.t20$geneSymbol), decreasing=TRUE)
#
# ELF1  ELK4 FOXO1 GABPA MEF2D 
#    1     1     1     1     1 
tbl.breaks.strong.t20  <- subset(tbl.breaks.t20, pctRef > 0.8)
sort(table(tbl.breaks.strong.t20$geneSymbol), decreasing=TRUE)
# GABPA   ELK4 NFATC3   ELF1  FOXO1  MEF2D  PRRX1   RBPJ   TBR1 
#     3      2      2      1      1      1      1      1      1 


print(load("ptk2b.eqtls.100.motifbreaks.RData"))
results.b100 <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.b100 <- tbl.breaks
tbl.summary.b100 <- tbl.summary
tbl.breaks.strong.b100  <- subset(tbl.breaks.b100, pctRef > 0.9)
sort(table(tbl.breaks.strong.b100$geneSymbol), decreasing=TRUE)
#   NFYC     ZFX   MEIS3   PRRX1    TCF3     MGA NEUROD1 NEUROD2    NFIA   NFKB1    NFYA     SP1     SP3    TBR1 
#      3       3       2       2       2       1       1       1       1       1       1       1       1       1 
tbl.breaks.strong.b100  <- subset(tbl.breaks.b100, pctRef > 0.8)
sort(table(tbl.breaks.strong.b100$geneSymbol), decreasing=TRUE)
#  MEIS3    NFIA   PRRX1    TCF3    IRF2    NFYC    RBPJ   GABPA     MGA    NFYA     SP1     SP3    TBR1     ZFX    FOSB NEUROD1 NEUROD2  NFATC3   SMAD2    ATF3   BACH1    E2F3    ELK4    MAFK  NFE2L2   NFKB1    RARA 
#      8       8       8       6       4       4       4       3       3       3       3       3       3       3       2       2       2       2       2       1       1       1       1       1       1       1       1 
#   RELA     SP4 
#      1       1 

rsids.b20 <- tail(tbl.eQTL$rsid, n=20)
tbl.summary.b20 <- tbl.summary.b20[, rsids.b20]
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.b20  <- subset(tbl.breaks.b100, pctRef > 0.9)
dim(tbl.breaks.b20)
dim(tbl.summary.b20)
dim(tbl.summary.t20)

goi <- c("GABPA", "PKNOX2", "MEF2C", "STAT4", "ELK4")

print(load("ptk2b.eqtls.top.200.motifbreaks.RData"))
results.t200 <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
dim(tbl.breaks)
tbl.breaks.t200 <- tbl.breaks
tbl.summary.t200 <- tbl.summary
tbl.breaks.strong.t200  <- subset(tbl.breaks.t200, pctRef > 0.9)
sort(table(tbl.breaks.strong.t200$geneSymbol), decreasing=TRUE)

#  ELF1  GABPA  MEIS3   RBPJ  STAT4    ZFX   ELK4   FOSB  FOXO1  MEF2C  MEF2D    MGA   MXI1 NFATC3   NFIA   NFYA   NFYC  SMAD2    SP1    SP3   TBR1 
#     2      2      2      2      2      2      1      1      1      1      1      1      1      1      1      1      1      1      1      1      1 
tbl.breaks.strong.t200  <- subset(tbl.breaks.t200, pctRef > 0.8)
sort(table(tbl.breaks.strong.t200$geneSymbol), decreasing=TRUE)
#  PRRX1    MXI1    TCF3   GABPA     MGA  NFATC3    RBPJ     SP1     SP3    ELF1    FOSB    HEY2   MEIS3   STAT4    ATF3    ELK4   FOXO1    MAFK   MEF2C NEUROD2    NFIA   SMAD2    EGR4   MEF2D NEUROD1  NFE2L2    TBR1 
#     10       9       9       6       6       6       6       6       6       5       5       5       5       5       4       4       4       4       4       4       4       4       3       3       3       3       3 
# ZBTB7B    IRF2   NFKB1    NFYA    TFEB    TP73     ZFX   BACH1 BHLHE41    CTCF    E2F3    NFYC   PPARG    RELA     SP4    TAF1 
#      3       2       2       2       2       2       2       1       1       1       1       1       1       1       1       1 




print(load("ptk2b.eqtls.random.200.06sep.1254.motifbreaks.RData"))
results.r200a <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.r200a <- tbl.breaks
tbl.summary.r200a <- tbl.summary
tbl.breaks.strong.r200a  <- subset(tbl.breaks.r200a, pctRef > 0.9)
sort(table(tbl.breaks.strong.r200a$geneSymbol), decreasing=TRUE)
#  FOXO1   MEIS3    TCF3     MGA  NFATC3    NFYC     SP3     ZFX    ELF1   GABPA    MAFK   MEF2D    MXI1 NEUROD2   PRRX1     SP1   STAT4 
#      4       3       3       2       2       2       2       2       1       1       1       1       1       1       1       1       1 

print(load("ptk2b.eqtls.random.200.06sep.1334.motifbreaks.RData"))
results.r200b <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.r200b <- tbl.breaks
tbl.summary.r200b <- tbl.summary
tbl.breaks.strong.r200b  <- subset(tbl.breaks.r200b, pctRef > 0.9)
sort(table(tbl.breaks.strong.r200b$geneSymbol), decreasing=TRUE)
#   TCF3     MGA    MXI1     ZFX  NFATC3    TBR1    E2F3    ELF1    HEY2   MEIS3 NEUROD2   PPARG   PRRX1   SMAD2     SP1   STAT4    TFEB 
#      4       3       3       3       2       2       1       1       1       1       1       1       1       1       1       1       1 


#    TCF3    MXI1    ELF1     MGA  NFATC3   PPARG     ZFX    E2F3    HEY2   MEIS3    TBR1    ELK4 NEUROD2   PRRX1   SMAD2     SP1   STAT4    TFEB 
#       6       5       4       3       3       3       3       2       2       2       2       1       1       1       1       1       1       1 

f <- "ptk2b.eqtls.top.40.06Sep.1450.motifbreaks.RData"
print(load(f))
results.t40a <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.t40a <- tbl.breaks
tbl.summary.t40a <- tbl.summary
tbl.breaks.strong.t40a  <- subset(tbl.breaks.t40a, pctRef > 0.9)
sort(table(tbl.breaks.strong.t40a$geneSymbol), decreasing=TRUE)
# ELF1 GABPA MEF2D   SP1 STAT4 
#    1     1     1     1     1 

f <- "ptk2b.eqtls.top.40.06Sep.1506.motifbreaks.RData"
print(load(f))
results.t40b <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.t40b <- tbl.breaks
tbl.summary.t40b <- tbl.summary
tbl.breaks.strong.t40b  <- subset(tbl.breaks.t40b, pctRef > 0.8)
sort(table(tbl.breaks.strong.t40b$geneSymbol), decreasing=TRUE)

f <- "ptk2b.eqtls.top.40.06Sep.1514.motifbreaks.RData"
print(load(f))
results.t40c <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.t40c <- tbl.breaks
tbl.summary.t40c <- tbl.summary
tbl.breaks.strong.t40c  <- subset(tbl.breaks.t40c, pctRef > 0.9)
sort(table(tbl.breaks.strong.t40c$geneSymbol), decreasing=TRUE)
# ELF1 GABPA MEIS3   MGA  MXI1 SMAD2 STAT4  TBR1 
#    2     1     1     1     1     1     1     1 
tbl.breaks.strong.t40c  <- subset(tbl.breaks.t40c, pctRef > 0.8)
sort(table(tbl.breaks.strong.t40c$geneSymbol), decreasing=TRUE)
#  GABPA    ELF1    ELK4     MGA   PRRX1    TBR1    TCF3 BHLHE41    HEY2    MAFK   MEF2C   MEIS3    MXI1  NFATC3    NFIA    NFYA   SMAD2     SP1     SP3   STAT4    TP73 
#      3       2       2       2       2       2       2       1       1       1       1       1       1       1       1       1       1       1       1       1       1 

f <- "ptk2b.eqtls.top.40.06Sep.1521.motifbreaks.RData"
print(load(f))
results.t40d <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.t40d <- tbl.breaks
tbl.summary.t40d <- tbl.summary
tbl.breaks.strong.t40d  <- subset(tbl.breaks.t40d, pctRef > 0.9)
sort(table(tbl.breaks.strong.t40d$geneSymbol), decreasing=TRUE)
# ELF1  ELK4 MEIS3  NFYA  NFYC 
#    1     1     1     1     1 
tbl.breaks.strong.t40d  <- subset(tbl.breaks.t40d, pctRef > 0.8)
sort(table(tbl.breaks.strong.t40d$geneSymbol), decreasing=TRUE)
#  PRRX1   MEIS3     MGA NEUROD2   NFKB1    EGR4    ELF1    ELK4    FOSB   GABPA    MAFK   MEF2C  NFATC3  NFE2L2    NFIA    NFYA    NFYC    RELA     SP4    TBR1    TCF3 
#      5       2       2       2       2       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 


# all eQTLs pval <= 0.01
f <- "ptk2b.eqtls.top.40.06Sep.1540.motifbreaks.RData"
print(load(f))
results.t45e <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.t45e <- tbl.breaks
tbl.summary.t45e <- tbl.summary
tbl.breaks.strong.t45e  <- subset(tbl.breaks.t45e, pctRef > 0.9)
sort(table(tbl.breaks.strong.t45e$geneSymbol), decreasing=TRUE)
# ELF1  ELK4  FOSB FOXO1 GABPA MEF2D  NFYA  NFYC STAT4 
#    2     1     1     1     1     1     1     1     1 
tbl.breaks.strong.t45e  <- subset(tbl.breaks.t45e, pctRef > 0.8)
sort(table(tbl.breaks.strong.t45e$geneSymbol), decreasing=TRUE)
#  GABPA    ELK4  NFATC3   PRRX1   STAT4    ELF1    FOSB NEUROD2    TCF3    EGR4   FOXO1    HEY2   MEF2C   MEF2D   MEIS3    MXI1 NEUROD1  NFE2L2    NFYA    NFYC    RBPJ     SP1     SP3     SRF    TBR1    TFEB 
#      4       3       3       3       3       2       2       2       2       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 


f <- "ptk2b.eqtls.middle.40.06Sep.1608.motifbreaks.RData"
print(load(f))
results.tmiddle.a <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.tmiddle.a <- tbl.breaks
tbl.summary.tmiddle.a <- tbl.summary
tbl.breaks.strong.tmiddle.a  <- subset(tbl.breaks.tmiddle.a, pctRef > 0.9)
sort(table(tbl.breaks.strong.tmiddle.a$geneSymbol), decreasing=TRUE)
#  FOXO1   MEIS3 NEUROD2  NFE2L2   PRRX1    RBPJ     SP3 
#      1       1       1       1       1       1       1 
tbl.breaks.strong.tmiddle.a  <- subset(tbl.breaks.tmiddle.a, pctRef > 0.8)
sort(table(tbl.breaks.strong.tmiddle.a$geneSymbol), decreasing=TRUE)
#   NFIA   MEIS3  NFE2L2    RBPJ    EGR3    ELF1    FOSB   FOXO1     MGA NEUROD2  NFATC3   NFKB1   PRRX1    RARA    RELA     SP1     SP3     SRF   STAT4    TBR1  ZBTB7B 
#      5       2       2       2       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1       1 

f <- "ptk2b.eqtls.middle.1001.06Sep.2334.motifbreaks.RData"
print(load(f))
results.middle1k <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.middle1k <- tbl.breaks
tbl.summary.middle1k <- tbl.summary
tbl.breaks.strong.middle1k  <- subset(tbl.breaks.middle1k, pctRef > 0.9)
sort(table(tbl.breaks.strong.middle1k$geneSymbol), decreasing=TRUE)
#   TCF3  NFATC3    NFYC   FOXO1   MEIS3   PRRX1   STAT4    RBPJ    ELF1    NFYA     ZFX    ELK4   GABPA     MGA  NFE2L2   NFKB1     SP1     SP3    ATF3    NFIA   SMAD2    TBR1   BACH1    E2F3    EGR3    EGR4    HEY2 
#     19      13      12      10       9       9       7       6       5       5       5       4       3       3       3       3       3       3       2       2       2       2       1       1       1       1       1 
#   MAFK   MEF2C   MEF2D    MXI1 NEUROD2    RARA    RELA     SRF    TAF1   TGIF2 
#      1       1       1       1       1       1       1       1       1       1 

tbl.breaks.strong.middle1k  <- subset(tbl.breaks.middle1k, pctRef > 0.8)
sort(table(tbl.breaks.strong.middle1k$geneSymbol), decreasing=TRUE)
#   TCF3    RBPJ    NFIA   PRRX1   MEIS3  NFATC3 NEUROD1   STAT4   FOXO1    MXI1    ELF1   GABPA     MGA    FOSB     SP1    ELK4 NEUROD2    TFEB   MEF2C    NFYA    NFYC     SP3    ATF3   SMAD2    HEY2    TBR1    E2F3 
#     62      61      48      47      42      37      33      29      27      27      26      25      23      22      20      19      19      19      18      18      18      17      16      15      14      14      13 
#    ZFX    MAFK  ZBTB7B  NFE2L2    TAF1 BHLHE41   MEF2D   NFKB1    RELA   BACH1    CTCF    IRF2    RARA     SRF    EGR4  PKNOX2   PPARG     SP4    EGR3    REST   TGIF2    MTF1    TP73 
#     13      10      10       9       9       8       8       8       8       7       5       5       5       5       3       3       3       3       2       2       2       1       1 


f <- "ptk2b.eqtls.top.3.07Sep.0630.motifbreaks.RData"
print(load(f))
results.top3 <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.breaks.top3 <- tbl.breaks
dim(tbl.breaks.top3)   # 54 27
tbl.summary.top3 <- tbl.summary
dim(tbl.summary.top3)
tbl.breaks.strong.top3  <- subset(tbl.breaks.top3, pctRef > 0.9)
sort(table(tbl.breaks.strong.top3$geneSymbol), decreasing=TRUE)
tbl.breaks.strong.top3  <- subset(tbl.breaks.top3, pctRef > 0.8)
sort(table(tbl.breaks.strong.top3$geneSymbol), decreasing=TRUE)


f <- "ptk2b.eqtls.top.14.07Sep.0640.motifbreaks.RData"
print(load(f))
results.top14 <- results
dim(tbl.breaks)  # 321 26
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
dups <- which(duplicated(tbl.breaks$sig))
printf("dups: %d", length(dups))
tbl.breaks <- tbl.breaks[-dups,]
tbl.breaks.top14 <- tbl.breaks
dim(tbl.breaks.top14)   # 238 27
tbl.summary.top14 <- tbl.summary
dim(tbl.summary.top14)
dim(subset(tbl.breaks.top14, pctRef > 0.9))
tbl.breaks.strong.top14  <- subset(tbl.breaks.top14, pctRef > 0.9)
sort(table(tbl.breaks.strong.top14$geneSymbol), decreasing=TRUE)
new.order <- order(tbl.breaks.top14$scoreDelta, decreasing=FALSE)
tbl.breaks.top14.ordered <- tbl.breaks.top14[new.order,]

tbl.breaks.strong.top14  <- subset(tbl.breaks.top14, pctRef > 0.8)
sort(table(tbl.breaks.strong.top14$geneSymbol), decreasing=TRUE)

subset(tbl.breaks.top14, pctRef > 0.7)[, coi]



print(load("ptk2b.eqtls.top.14.08Sep.0947.motifbreaks.RData"))
results.t14 <- results
tbl.breaks$signature <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$signature))),]
tbl.breaks.t14 <- tbl.breaks
tbl.summary.t14 <- tbl.summary
tbl.breaks.strong.t14  <- subset(tbl.breaks.t14, pctRef > 0.9)
sort(table(tbl.breaks.strong.t14$geneSymbol), decreasing=TRUE)

tbl.breaks.strong.t14  <- subset(tbl.breaks.t14, pctRef > 0.8)
sort(table(tbl.breaks.strong.t14$geneSymbol), decreasing=TRUE)


new.order <- order(abs(tbl.breaks.t14$scoreDelta))

head(tbl.breaks.t14[new.order,], n=10)

# just two tfsPKNOX2 and GABPA

coi <- c("geneSymbol", "SNP_id", "pctRef", "pctAlt", "seqMatch", "pctDelta")
print(load("ptk2b.eqtls.top.14.08Sep.1017.motifbreaks.RData"))
results.t14 <- results
tbl.strong <- subset(as.data.frame(results, row.names=NULL), effect=="strong")
tbl.strong$pctDelta <- tbl.strong$pctAlt-tbl.strong$pctRef

tbl.strong[, coi]

tbl.breaks$signature <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$signature))),]
tbl.breaks.t14 <- tbl.breaks
tbl.summary.t14 <- tbl.summary
tbl.breaks.strong.t14  <- subset(tbl.breaks.t14, pctRef > 0.9)
sort(table(tbl.breaks.strong.t14$geneSymbol), decreasing=TRUE)

tbl.breaks.strong.t14  <- subset(tbl.breaks.t14, pctRef > 0.8)
sort(table(tbl.breaks.strong.t14$geneSymbol), decreasing=TRUE)


new.order <- order(abs(tbl.breaks.t14$scoreDelta))

head(tbl.breaks.t14[new.order,], n=10)


