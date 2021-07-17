exploreALFA <- function(rsid="rs372089992")
{
   library(httr)
   library(jsonlite)

   uri <- sprintf("https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/%s/frequency",
                  sub("rs", "", rsid))

   response <- GET(uri)
   suppressMessages(x <- fromJSON(httr::content(response, as="text"))$results)
   ref <- x[[1]]$ref
   counts.list <- x[[1]]$counts

   population.codes <- names(counts.list[[1]][[1]])
   tbls <- list()
   for(pop in population.codes){
       pop.record <- counts.list[[1]][[1]][[pop]]
       tbl <- data.frame(population=pop, ref=ref, rsid=rsid, stringsAsFactors=FALSE)
       alleles <- names(pop.record)
       counts <- as.integer(pop.record)
       tbl[1, alleles] <- counts
       tbls[[pop]] <- tbl
       }
   tbl.alfa <- do.call(rbind, tbls)
   rownames(tbl.alfa) <- NULL

   f <- system.file(package="EndophenotypeExplorer", "extdata", "gwas", "populationCodes.tsv")
   stopifnot(file.exists(f))
   tbl.popCodes <- read.table(f, header=TRUE, as.is=TRUE, sep="\t")

   tbl.alfa <- merge(tbl.alfa, tbl.popCodes, by.x="population", by.y="biosampleID")
   allele.colnames <- setdiff(colnames(tbl.alfa), c("population", "ref", "name", "rsid"))
   tbl.alfa <- tbl.alfa[, c("rsid", "ref", "name", allele.colnames)]
   deleter <- grep("population", colnames(tbl.alfa))
   if(length(deleter) == 1)
       tbl.alfa <- tbl.alfa[, -(deleter)]
   coi <- c("rsid", "ref", "name", allele.colnames)
   tbl.alfa <- tbl.alfa[, coi]
   colnames(tbl.alfa)[grep("name", colnames(tbl.alfa))] <- "population"
   tbl.alfa$total <- rowSums(tbl.alfa[, allele.colnames])

   for(allele in allele.colnames){
      allele.freq <- 100 * tbl.alfa[[allele]]/tbl.alfa$total
      new.colname <- sprintf("%s.freq", allele)
      tbl.alfa[[new.colname]] <- allele.freq
      } # for allele

   tbl.alfa

} # exploreALFA
#----------------------------------------------------------------------------------------------------
test_alfa <- function()
{

    message(sprintf("--- test_alfa"))

        # these are all associated with the BIN1 GWAS landmark snp according to
        # Posthuma Genome-wide meta-analysis identifies new loci and functional pathways
        # influencing  Alzheimer's disease risk
        #   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6836675/

    rsids <- c("rs114469753",
               "rs11674295",
               "rs11674920",
               "rs11675014",
               "rs11680116",
               "rs11680117",
               "rs12472671",
               "rs12475915",
               "rs12477427",
               "rs12617835",
               "rs12618593",
               "rs13019485",
               "rs202224448",
               "rs2118506",
               "rs2118508",
               "rs34558561",
               "rs34997561",
               "rs35103166",
               "rs35114168",
               "rs35155471",
               "rs35400498",
               "rs35486865",
               "rs35783077",
               "rs372089992",
               "rs3943703",
               "rs4663105",
               "rs60447541",
               "rs66837244",
               "rs6710467",
               "rs7561528")

    set.seed(17)
       # rs184384746 is rare in Europeans, worth a look
    rsids <- c("rs184384746", sort(sample(rsids, size=5)))
    tbls <- list()

    for(rsid in rsids){
       printf("--- %s", rsid)
       tbl <- exploreAlfa(rsid)
       tbls[[rsid]] <- tbl
       Sys.sleep(1)
       }

   tbl <- do.call(rbind.fill, tbls)
   checkEquals(dim(tbl), c(60, 12))

      # an idiosyncratic test
   tbl.af.rare <- subset(tbl, population=="African American" & (A.freq < 1 | G.freq <1 | C.freq <1 | T.freq < 1))
   checkEquals(tbl.af.rare$rsid, "rs114469753")
      #          rsid ref       population C    T total  C.freq   T.freq  A  G A.freq G.freq
      #   rs114469753   T African American 5 2851  2856 0.17507 99.82493 NA NA     NA     NA


} # test_alfa
#----------------------------------------------------------------------------------------------------
