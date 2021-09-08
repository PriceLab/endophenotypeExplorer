library(EndophenotypeExplorer)

etx <- EndophenotypeExplorer$new("PTK2B", "hg38");
tbl.eQTL <- etx$getEQTLsForGene()
tbl.eQTL  <- subset(tbl.eQTL, study=="ampad-mayo" & tissue=="tcx")
new.order <- order(tbl.eQTL$pvalue, decreasing=FALSE)
tbl.eQTL <- tbl.eQTL[new.order,]

dim(tbl.eQTL) # 5590 10
head(tbl.eQTL)

snpoi <- head(tbl.eQTL$rsid, n=100)
goi <- c("GABPA", "PKNOX2", "MEF2C", "STAT4", "ELK4")
coi <- c("geneSymbol", "SNP_id", "pctRef", "pctDelta", "seqMatch")

f <- "ptk2b.eqtls.top.200.motifbreaks.RData"
print(load(f))
results.t200 <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.best.09 <- subset(tbl.breaks, geneSymbol %in% goi & effect=="strong" & SNP_id %in% snpoi & pctRef > 0.9)[, coi]
tbl.best.09
#      geneSymbol     SNP_id    pctRef    pctAlt              seqMatch             pctDelta
# 17        GABPA  rs1000987 0.9011795 0.7454862       atagcacaGgaagatgaat       -0.15569330
# 2964       ELK4 rs73554189 0.9591010 0.8580483      ggcgtgaacCcggaaggcgga      -0.10105277
# 2967      GABPA rs73554189 0.9358740 0.8407787   aatggcgtgaacCcggaaggcggagct   -0.09509528
# 3047      STAT4 rs35491396 0.9672126 0.7818012      gtttcctccAaagttggccga      -0.18541141
# 4484      STAT4  rs7840424 0.9276379 0.7520905      agaatttgaGgaatcagtaaa      -0.17554734
match(unique(tbl.best.09$SNP_id), tbl.eQTL$rsid) # 16 24 67 44
tbl.eQTL[match(unique(tbl.best.09$SNP_id), tbl.eQTL$rsid),]
#       chrom     hg19     hg38       rsid      pvalue            ensg genesymbol      study tissue   assay
# 15643  chr8 26579782 26722265  rs1000987 0.005483287 ENSG00000120899      PTK2B ampad-mayo    tcx unknown
# 3739   chr8 26558310 26700793 rs73554189 0.006188634 ENSG00000120899      PTK2B ampad-mayo    tcx unknown
# 1843   chr8 26472064 26614548 rs35491396 0.024536124 ENSG00000120899      PTK2B ampad-mayo    tcx unknown
# 4215   chr8 27082640 27225123  rs7840424 0.009057078 ENSG00000120899      PTK2B ampad-mayo    tcx unknown

tbl.best.08 <- subset(tbl.breaks, geneSymbol %in% goi & effect=="strong" & SNP_id %in% snpoi & pctRef > 0.8)[, coi]
tbl.best.08
#      geneSymbol     SNP_id    pctRef    pctDelta               seqMatch
# 13         ELK4  rs1000987 0.8080655 -0.15248831             tagcacaGgaagatgaa              
# 17        GABPA  rs1000987 0.9011795 -0.15569330            atagcacaGgaagatgaat             
# 373       GABPA  rs6987195 0.8039934 -0.17064474        ccaagaaggcaaCgaaagtgacctctg         
# 459       MEF2C  rs2046222 0.8453355  0.06017737       atttaggtcaggtActatttctgtcattt        
# 916       MEF2C  rs7005443 0.8175658 -0.11786215         aatccagataaAtatttagaagctt          
# 930       STAT4  rs7005443 0.8533630 -0.18345971           tccagataaAtatttagaagc            
# 1869      STAT4 rs10503794 0.8836616 -0.16993282           ttgcttgctTctgtggatcgc            
# 2770      GABPA rs73554163 0.8057356 -0.12128561            tgccctctTcggggatgaa             
# 2964       ELK4 rs73554189 0.9591010 -0.10105277           ggcgtgaacCcggaaggcgga            
# 2967      GABPA rs73554189 0.9358740 -0.09509528        aatggcgtgaacCcggaaggcggagct         
# 3047      STAT4 rs35491396 0.9672126 -0.18541141           gtttcctccAaagttggccga            
# 4484      STAT4  rs7840424 0.9276379 -0.17554734           agaatttgaGgaatcagtaaa            
# 5275       ELK4 rs13281659 0.8015691 -0.07126471             gcttgaaTccgggaggt              
# 6707       ELK4 rs60565850 0.8973199 -0.17269081          acaagagattCcggtactaaaca           
# 6711      GABPA rs60565850 0.8424593 -0.17064474        gcacaagagattCcggtactaaacatt         
match(unique(tbl.best.08$SNP_id), tbl.eQTL$rsid)
#  [1] 16 85 98 22 34 19 24 67 44 60 20
length(unique(tbl.breaks$SNP_id)) # [1] 200
tbl.eQTL[match(unique(tbl.best.08$SNP_id), tbl.eQTL$rsid),]
#       chrom     hg19     hg38       rsid      pvalue
# 15643  chr8 26579782 26722265  rs1000987 0.005483287
# 3233   chr8 26524958 26667441  rs6987195 0.029883003
# 917    chr8 26212413 26354897  rs2046222 0.035056138
# 3395   chr8 26556825 26699308  rs7005443 0.006188634
# 15991  chr8 26556765 26699248 rs10503794 0.006188634
# 3738   chr8 26554960 26697443 rs73554163 0.006188589
# 3739   chr8 26558310 26700793 rs73554189 0.006188634
# 1843   chr8 26472064 26614548 rs35491396 0.024536124
# 4215   chr8 27082640 27225123  rs7840424 0.009057078
# 233    chr8 26526449 26668932 rs13281659 0.019897566
# 2858   chr8 26557329 26699812 rs60565850 0.006188634

f <- "ptk2b.eqtls.bottom.200.08Sep.1212.motifbreaks.RData"
print(load(f))
results.b200 <- results
tbl.breaks$sig <- with(tbl.breaks, sprintf("%s:%s", geneSymbol, SNP_id))
tbl.breaks <- tbl.breaks[-(which(duplicated(tbl.breaks$sig))),]
tbl.bottom.200.09 <- subset(tbl.breaks, geneSymbol %in% goi & effect=="strong" & pctRef > 0.9)[, coi]
tbl.bottom.200.09
# nothing

tbl.bottom.200.08 <- subset(tbl.breaks, geneSymbol %in% goi & effect=="strong" & pctRef > 0.8)[, coi]
tbl.bottom.200.08
#     geneSymbol    SNP_id    pctRef   pctDelta                    seqMatch
# 123      GABPA rs2292976 0.8392180  0.1016270 tcctccccttccggaagagCccttat 
# 287      GABPA rs3735758 0.8365278 -0.1669636 tggtagaggggaaggggctCatttga 

match(unique(tbl.bottom.200.08$SNP_id), tbl.eQTL$rsid) # [1] 5494 5492

tbl.eQTL[match(unique(tbl.bottom.200.08$SNP_id)
#    chrom     hg19     hg38      rsid    pvalue 
# 44  chr8 27318544 27461027 rs2292976 0.9847228 
# 26  chr8 27310546 27453029 rs3735758 0.9842101 

match(unique(tbl.bottom.200.08$SNP_id), tbl.eQTL$rsid)
[1] 5494 5492
