# rs4575098
# allegedly associated with  PPOX,APOA2,ADAMTS4,B4GALT3,USP21,NDUFS2,TOMM40
# 1:161185602 (GRCh38)
# 1:161155392 (GRCh37)
#----------------------------------------------------------------------------------------------------
library(EndophenotypeExplorer)
library(plyr)
library(ghdb)
#----------------------------------------------------------------------------------------------------
ghdb <- GeneHancerDB()
tag.snp <- "rs4575098"
tbl.williams <- get(load(system.file(package="EndophenotypeExplorer", "extdata", "gwas",
                                     "williams-natureNeuroscience2020.RData")))
dim(tbl.williams)
rsids.williams <- unique(grep("^rs", unlist(strsplit(tbl.williams$rsid, ",")), value=TRUE))
tag.snp %in% rsids.williams
tbl.posthuma <- get(load(system.file(package="EndophenotypeExplorer", "extdata", "gwas",
                                     "tbl.posthuma-38-geneAssociations-curated-3828x12.RData")))
dim(tbl.posthuma)
tag.snp %in% tbl.posthuma$rsid

subset(tbl.williams, rsid==tag.snp)
#    locusOrGene      rsid
# 54     ADAMTS4 rs4575098

igv <- start.igv("ADAMTS4", "hg19")
zoomOut(igv)
