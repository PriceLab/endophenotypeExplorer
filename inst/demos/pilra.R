# request from cory, email (23 jul 2021) "PILRA data query"
# From the VCF files, we can identify all AMP-AD samples with the rs1859788 variant (Paul).a

# PILRA is a protective AD GWAS variant that we believe is tied to activation of microglia
# cells. There is a protective mutation in G78R, or rs1859788, that occurs at an allele frequency of
# ~0.29 in the population.

# From the VCF files, we can identify all AMP-AD samples with the rs1859788 variant (Paul).

# For the ROSMAP data, we have CERAD scores, but I think data also exists for things like A-beta
# levels and p-Tau. Tain: I think you may have asked Vilas about this. We'd love to get more
# phenotypic data for ROSMAP than we currently have.

# From the ROSMAP data, we have the percentages of microglial cells for the subpopulations (Andrew).

# There are two questions I'd like to ask:

# 1) Is there a difference in microglial subpopulation percentages in the rs1859788 containing
# samples? (Paul, Andrew)

# 2) Is there a difference in phenotypic endpoints (CERAD, or whatever additional phenotypic data we
# have for ROSMAP) in rs1859788 containing samples.

# For both these questions, I think it makes sense to do both AD vs control, and another analysis
# with both AD and control.
#
#----------------------------------------------------------------------------------------------------
library(EndophenotypeExplorer)
tag.snp <- "rs1859788"
etx <- EndophenotypeExplorer$new("PILRA", "hg19")
etx$getAggregatedAlleleFrequencies(tag.snp)

#         rsid ref       population C     G     A total C.freq   G.freq   A.freq min.freq
# 1  rs1859788   A         European 0 44875 20233 65108      0 68.92394 31.07606 31.07606
# 2  rs1859788   A   African Others 0    52     6    58      0 89.65517 10.34483 10.34483
# 3  rs1859788   A       East Asian 0    68   122   190      0 35.78947 64.21053 35.78947
# 4  rs1859788   A African American 0  1653   381  2034      0 81.26844 18.73156 18.73156
# 5  rs1859788   A Latin American 1 0   345   155   500      0 69.00000 31.00000 31.00000
# 6  rs1859788   A Latin American 2 0   625   825  1450      0 43.10345 56.89655 43.10345
# 7  rs1859788   A      Other Asian 0    41    55    96      0 42.70833 57.29167 42.70833
# 8  rs1859788   A      South Asian 0     9    27    36      0 25.00000 75.00000 25.00000
# 9  rs1859788   A          African 0  1705   387  2092      0 81.50096 18.49904 18.49904
# 10 rs1859788   A            Asian 0   109   177   286      0 38.11189 61.88811 38.11189
# 11 rs1859788   A            Total 0 53386 24396 77782      0 68.63542 31.36458 31.36458
# 12 rs1859788   A            Other 0  5718  2592  8310      0 68.80866 31.19134 31.19134
etx$rsidToLoc(tag.snp)

#   chrom     hg19      hg38      rsid
# 1     7 99971834 100374211 rs1859788

loc <- 99971834
mtx.geno <- etx$getGenoMatrix("7", loc-1, loc+1)
dim(mtx.geno)
table(mtx.geno[1,])
# 0/0 0/1 1/1
# 185 800 909
