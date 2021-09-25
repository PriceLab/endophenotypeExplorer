library(EndophenotypeExplorer)
library(RUnit)
#----------------------------------------------------------------------------------------------------
etx <- EndophenotypeExplorer$new(NA, NA, initialize.snpLocs=TRUE)
tbl.gwas <- etx$getGWASTables()$gwas.38
#----------------------------------------------------------------------------------------------------
lookupGenotypesWriteTable <- function()
{
   geno.matrices <- list()
   for(r in seq_len(nrow(tbl.gwas))){
      gene <- tbl.gwas$geneSymbol[r]
      rsid <- tbl.gwas$leadVariant[r]
      etx$setTargetGene(gene, "hg19")
      printf("--- getting genoMatrix for %s", rsid)
      mtx.geno <- etx$getGenoMatrixByRSID(rsid)
      geno.matrices[[rsid]] <- mtx.geno
      }

   rsids <- tbl.gwas$leadVariant
   length(rsids) # 38
   mtx.genos <- do.call(rbind, geno.matrices)
   dim(mtx.genos)
   tbl.genos <- as.data.frame(mtx.genos, check.names=FALSE)
   rownames(tbl.genos) <- rsids
   dim(tbl.genos)
   tbl.out <- cbind(tbl.gwas, tbl.genos)
   grep("leadVariant", colnames(tbl.out))
   tbl.out <- tbl.out[, -2]
   write.table(tbl.out, file="adGwas38-ampad-genotypes.tsv", sep="\t", quote=FALSE,
               row.names=TRUE, col.names=TRUE)

} # lookupGenotypesWriteTable
#----------------------------------------------------------------------------------------------------
test_table <- function()
{
    tbl <- read.table("adGwas38-ampad-genotypes.tsv", header=TRUE,
                      check.names=FALSE, sep="\t", as.is=TRUE)
    checkEquals(dim(tbl), c(38, 1906))

      # choose some snps at random, make sure that the genotypes agree with
      # those returned by a fresh call

    rows <- sort(sample(seq_len(nrow(tbl)), size=3))
    for(row in rows){
      targetGene <- tbl.gwas$geneSymbol[row]
      etx$setTargetGene(targetGene, "hg19")
      mtx.geno <- etx$getGenoMatrixByRSID(tbl.gwas$leadVariant[row])
      checkEquals(as.character(mtx.geno[1,]), as.character(tbl[row, 13:1906]))
      }

} # test_table
#----------------------------------------------------------------------------------------------------



