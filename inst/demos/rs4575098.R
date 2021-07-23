# rs4575098,
# 0.98 LD with 	rs11585858
#   chrom      hg19      hg38       rsid
# 1     1 161156033 161186243 rs11585858
# allegedly associated with  PPOX,APOA2,ADAMTS4,B4GALT3,USP21,NDUFS2,TOMM40
# 1:161185602 (GRCh38)
# 1:161155392 (GRCh37)
#
# not reported by posthuma
# from williams:
#   ADAMTS4 rs4575098; ref. 31
#   Combined UK Biobank, ADSP, IGAP, PGC–ALZ and deCODE
#   Extracellular matrix metalloproteinase (aggrecanase-1)
#
# Jansen, I. E. et al. Genome-wide meta-analysis identifies new loci and
# functional pathways influencing Alzheimer’s disease risk. Nat. Genet. 51,
# 404–413 (2019).  29 loci in total, implicating 215 potential causative genes
# here is our tag snp
#
#          region                     case-control status(phase 1)        AD-by-proxy (phase 2)                 overall (phase3)
#   locus    chr  loc         gene          snp     p                             SNP    p               A1 A2 MAF    Z       p
#    1       chr1 161155392  ADAMTS4      rs457098  1.57e-04                  rs457098   6.88e-08         A  G 0.240 6.36  2.05E-10
#
# Here, we performed a large genome-wide association study of clinically diagnosed AD and
# AD-by-proxy (71,880 cases, 383,378 controls). AD-by- proxy, based on parental diagnoses, showed
# strong genetic correlation with AD (rg=0.81). Meta- analysis identified 29 risk loci, implicating
# 215 potential causative genes. Associated genes are strongly expressed in immune-related tissues
# and cell types (spleen, liver and microglia). Gene-set analyses indicate biological mechanisms
# involved in lipid-related processes and degradation of amyloid precursor proteins. We show strong
# genetic correlations with multiple health-related outcomes, and Mendelian randomisation results
# suggest a protective effect of cognitive ability on AD risk. These results are a step forward in
# identifying the genetic factors that contribute to AD risk and add novel insights into the
# neurobiology of AD.
# new table here, once published:
#
# Bellenguez C, Küçükali F, Jansen I, et al. New insights on the genetic etiology of Alzheimer’s
# and related dementia. medRxiv; 2020. DOI: 10.1101/2020.10.01.20200659.

#----------------------------------------------------------------------------------------------------
library(EndophenotypeExplorer)
library(plyr)
library(ghdb)
#----------------------------------------------------------------------------------------------------
ghdb <- GeneHancerDB()
tag.snp <- "rs4575098"
ld.snp  <- "rs11585858"
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

etx <- EndophenotypeExplorer$new("ADAMTS4", "hg19")

tbl.loc <- etx$rsidToLoc(tag.snp)
   #   chrom      hg19      hg38      rsid
   #       1 161155392 161185602 rs4575098

tbl.track <- with(tbl.loc, data.frame(chrom=chrom, start=hg19-1, end=hg19, stringsAsFactors=FALSE))
track <- DataFrameAnnotationTrack(tag.snp, tbl.track, color="red", trackHeight=25)
displayTrack(igv, track)

#   chrom      hg19      hg38       rsid
# 1     1 161156033 161186243 rs11585858
tbl.track <- data.frame(chrom="1", start=161156033-1, end=161156033, stringsAsFactors=FALSE)
track <- DataFrameAnnotationTrack(ld.snp, tbl.track, color="red", trackHeight=25)
displayTrack(igv, track)



tbl.all <- retrieveEnhancersFromDatabase(ghdb, target.gene, tissues="all")
dim(tbl.all) # 71 16
tbl.hg19 <- to.hg19(ghdb, tbl.all)
dim(tbl.hg19)
tbl.track <- tbl.hg19[, c("chrom", "start", "end", "combinedscore")]
colnames(tbl.track)[4] <- "score"
tbl.track$score <- asinh(tbl.track$score)
track <- DataFrameQuantitativeTrack("ADAMST4-gh", tbl.track, color="blue", autoscale=TRUE)
#displayTrack(igv, track)


# hg38:  chr1, 161181277 161188841
tbl.ghByRegion <- queryByRegion(ghdb, chrom="1", start=161181277, end=161188841)
dim(tbl.ghByRegion)
tbl.ghByRegion

genes.other <- unique(tbl.ghByRegion$gene)
for(gene in genes.other){
   tbl.gh <- retrieveEnhancersFromDatabase(ghdb, gene, tissues="all")
   tbl.hg19 <- to.hg19(ghdb, tbl.gh)
   tbl.track <- tbl.hg19[, c("chrom", "start", "end", "combinedscore")]
   colnames(tbl.track)[4] <- "score"
   tbl.track$score <- asinh(tbl.track$score)
   trackName <- sprintf("%s-gh", gene)
   track <- DataFrameQuantitativeTrack(trackName, tbl.track, color="random", autoscale=FALSE,
                                       min=0, max=10)
   displayTrack(igv, track)
   }


#----------------------------------------------------------------------------------------------------
explore.variants <- function()
{
   showGenomicRegion(igv, "chr1:161,155,310-161,156,114")  # tag & ld snp, with one inbetween
     # maybe adjust manually, then:
   roi.igv <- getGenomicRegion(igv)
   roi.igv$chrom <- sub("chr", "", roi.igv$chrom)
   #roi <- with(roi.igv, GRanges(seqnames=chrom, IRanges(start=start, end=end)))
   #url <- "https://igv-data.systemsbiology.net/static/ampad/NIA-1898/chr1.vcf.gz"

   mtx.geno <- with(roi.igv, etx$getGenoMatrix(chrom, start, end))
   dim(mtx.geno)
   rsid <- etx$locsToRSID(rownames(mtx.geno), "hg19")
   rownames(mtx.geno) <- rsid
   rownames(mtx.geno)

   non.wt <- apply(mtx.geno, 1, function(row) length(which(row != "0/0")))
   tbls <- list()
   rsids.legit <- grep("^rs", names(non.wt), value=TRUE)
   for(rsid in rsids.legit){
      printf("--- rsid: %s", rsid)
      tbls[[rsid]] <- etx$getAggregatedAlleleFrequencies(rsid)
      }
   tbl <- do.call(rbind.fill, tbls)
   tbl.euro <- subset(tbl, population=="European")

   tbl.euro
       #          rsid ref population     A      G  total      A.freq    G.freq     min.freq     C    C.freq    T   T.freq
      ##13   rs4575098   G   European 40918 129876 170794 23.95751607  76.04248  23.95751607    NA        NA   NA       NA
      ##73  rs11585858   C   European 27919      0 121662 22.94800349   0.00000  22.94800349 93743  77.05200   NA       NA
      ##49  rs60677420   C   European    NA     NA  24808          NA        NA   7.20735247 23020  92.79265 1788 7.207352

      # 1  rs774377891  G   European    13  14273  14286  0.09099818  99.90900   0.09099818    NA        NA   NA       NA
      #25 rs140964134   C   European     2      0  14266  0.01401935   0.00000   0.01401935 14264  99.98598   NA       NA
      #37 rs181507083   G   European     0   9616   9616  0.00000000 100.00000 100.00000000    NA        NA    0 0.000000
      #61 rs764488194   C   European    NA     NA   9690          NA        NA 100.00000000  9690 100.00000    0 0.000000

    table(mtx.geno["rs4575098",])
      #  0/0  0/1  1/1      40%
      # 1121  682   91
    table(mtx.geno["rs11585858",])
      #  0/0  0/1  1/1      40%
      # 1120  681   93
    table(mtx.geno["rs60677420", ])
      #  0/0  0/1  1/1
      # 1600  280   14      15%

    samples.mutant.rs4575098 <- colnames(mtx.geno)[which(mtx.geno["rs4575098",] != "0/0")]
    samples.mutant.rs11585858 <- colnames(mtx.geno)[which(mtx.geno["rs11585858",] != "0/0")]
    samples.mutant.rs60677420 <- colnames(mtx.geno)[which(mtx.geno["rs60677420",] != "0/0")]
    length(samples.mutant.rs4575098)    # 773
    length(samples.mutant.rs11585858)   # 774
    length(samples.mutant.rs60677420)   # 294

    length(intersect(samples.mutant.rs4575098, samples.mutant.rs11585858))   # 768
    length(intersect(samples.mutant.rs60677420, samples.mutant.rs4575098))   # 55

     # up next: are these all really european?  all path aging?  all with the variant
     # shown above, at 23% euro frequency? map to patient, examine ethnicity
     # more flexible rownames look up for mtx.geno.  they originally had the substitution
     # is that preserved in the rsid, or covered over when there are 3 or 4 alleles?

     # a small case to start, the 14 homozygous samples for rs60677420
     soi <- "rs4575098"
     soi <- "rs11585858"
     soi <- "rs60677420"
     vcf.ids <- names(which(mtx.geno[soi,] != "0/0"))
     length(vcf.ids)
     tbl.map <- etx$getIdMap()
     stopifnot(all(vcf.ids %in% tbl.map$vcf))
     tbls.clinical <- lapply(vcf.ids, function(id) etx$vcfSampleID.to.clinicalTable(id))
     tbl.clinical <- do.call(rbind, tbls.clinical)
     mean(tbl.clinical$cogdx, na.rm=TRUE) # [1] 2.627957
     mean(tbl.clinical$braak, na.rm=TRUE) # [1] 3.779539
     mean(tbl.clinical$cerad, na.rm=TRUE) # [1] 2.207705

       #---------------------------------
       # stats for all samples
       #---------------------------------
     tbls.clinical <- lapply(colnames(mtx.geno), function(id) etx$vcfSampleID.to.clinicalTable(id))
     tbl.clinical <- do.call(rbind, tbls.clinical)
     mean(tbl.clinical$cogdx, na.rm=TRUE) # [1] 2.653913
     mean(tbl.clinical$braak, na.rm=TRUE) # [1] 3.65516
     mean(tbl.clinical$cerad, na.rm=TRUE) # [1] 2.204013


       #--------------------------------------------------
       # get mtx.ad
       #--------------------------------------------------
     samples.ad <- subset(tbl.clinical, cogdx >= 4)$sampleID
     length(samples.ad)  # 498
     snps.oi <- c("rs4575098", "rs11585858", "rs60677420")
     mtx.ad <- mtx.geno[snps.oi, samples.ad]
     dim(mtx.ad)   # 3 498
     table(mtx.ad[1,])
     table(mtx.ad[2,])
     table(mtx.ad[3,])

       #----------------------------------------
       # pvalue of 098 in
       #----------------------------------------

   studyPop.rs4575098 <- rep(0, ncol(mtx.ad))
   studyPop.rs4575098[which(mtx.ad[1,] != "0/0")] <- 1
   table(studyPop.rs4575098)

   euroPop.all <- rep(0, 170794)
   euroPop.all[sample(seq_len(length(euroPop.all)), 40918)] <- 1
   table(euroPop.all)
   t.test(studyPop.rs4575098, euroPop.all)$p.value  # [1] 4.178421e-12

} # explore.variants
#----------------------------------------------------------------------------------------------------
# cribbed from https://github.com/KatrionaGoldmann/omicAnnotations/blob/main/R/gtex_eqtl.R
# see also: https://smart-api.info/ui/8b3389ee427f89a358e157319a9db534#/expression/geneExpression
# gene id mapping: https://cran.r-project.org/web/packages/grex/vignettes/grex.html
# Why does GTEx use TPM units rather than RPKM?
# We no longer provide expression quantifications in RPKM. TPM is a better unit for comparison of
# RNA-seq samples. You can convert RPKM as follows (for each sample/column):
# TPM = RPKM / sum(RPKM) * 1e6. For additional information,
# please see: https://academic.oup.com/bioinformatics/article/26/4/493/243395
gtex.demo <- function()
{
   library(httr)
   library(grex); data(gtexv7)
   url <-  paste0("https://gtexportal.org/rest/v1/reference/gene?geneId=", 10,
                  "&gencodeVersion=v26&genomeBuild=GRCh38%2Fhg38&pageSize=250",
                  "&format=json")
   x <- GET(url)
   content(x)

   url.0 <- "https:///v1/expression/geneExpression?"
   url.1 <- "datasetId=gtex_v8&gencodeId=ENSG00000065613.9"
   url.2 <- "&tissueSiteDetailId=Brain_Cortex&format=json"

   url <- sprintf("%s%s%s", url.0, url.1, url.2)
   x <- GET(url)
   content(x)
   url <- https://gtexportal.org/rest/v1/expression/geneExpression?datasetId=gtex_v8&gencodeId=ENSG00000065613.9&tissueSiteDetailId=Brain_Cortex&format=json

   tissues <- c('Adipose_Subcutaneous', 'Adipose_Visceral_Omentum',
                'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary',
                'Artery_Tibial', 'Bladder', 'Brain_Amygdala',
                'Brain_Anterior_cingulate_cortex_BA24',
                'Brain_Caudate_basal_ganglia',
                'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum',
                'Brain_Cortex', 'Brain_Frontal_Cortex_BA9',
                'Brain_Hippocampus', 'Brain_Hypothalamus',
                'Brain_Nucleus_accumbens_basal_ganglia',
                'Brain_Putamen_basal_ganglia',
                'Brain_Spinal_cord_cervical_c-1',
                'Brain_Substantia_nigra', 'Breast_Mammary_Tissue',
                'Cells_EBV-transformed_lymphocytes',
                'Cells_Cultured_fibroblasts',
                'Cells_Transformed_fibroblasts', 'Cervix_Ectocervix',
                'Cervix_Endocervix', 'Colon_Sigmoid',
                'Colon_Transverse',
                'Esophagus_Gastroesophageal_Junction',
                'Esophagus_Mucosa', 'Esophagus_Muscularis',
                'Fallopian_Tube', 'Heart_Atrial_Appendage',
                'Heart_Left_Ventricle', 'Kidney_Cortex',
                'Kidney_Medulla', 'Liver', 'Lung',
                'Minor_Salivary_Gland', 'Muscle_Skeletal',
                'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary',
                'Prostate', 'Skin_Not_Sun_Exposed_Suprapubic',
                'Skin_Sun_Exposed_Lower_leg',
                'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach',
                'Testis', 'Thyroid', 'Uterus', 'Vagina',
                'Whole_Blood')

   gencode.gr <- import("gencode.v26.GRCh38.genes.gtf")
   tbl.gencode <- as.data.frame(gencode.gr)

   brain.tissues <- grep("Brain", tissues, value=TRUE)
   brain.tissue.string <- paste(brain.tissues, collapse=",")
   gata2 <- "ENSG00000179348.11"
   fcer1g <- "ENSG00000158869.10"

   geneSymbol <- "ADAMTS4"
   geneSymbol <- "FCER1G"
   gencode.id <- subset(tbl.gencode, type=="gene" & transcript_name==geneSymbol)$gene_id

   url <- "https://gtexportal.org/rest/v1/expression/geneExpression?datasetId=gtex_v7&gencodeId=ENSG00000135100.13&tissueSiteDetailId=Pancreas&format=json"
   GET(url)

   url.0 <- "https://gtexportal.org/rest/v1/expression/geneExpression"
   url.0 <- "https://gtexportal.org/rest/v1/expression/medianGeneExpression"
   url.1 <- "datasetId=gtex_v7"
   url.1 <- "datasetId=gtex_v8"
   url.2 <- sprintf("gencodeId=%s", gata2)
   url.2 <- sprintf("gencodeId=%s", fcer1g)
   url.2 <- sprintf("gencodeId=%s", gencode.id)
   #url.2 <- sprintf("gencodeId=%s", "ENSG00000135100")
   #url.3 <- "tissueSiteDetailId=Pancreas,Brain_Cerebellum"
   url.3 <- sprintf("tissueSiteDetailId=%s", brain.tissue.string)
   url.4 <- "format=json"
   url <- sprintf("%s?%s&%s&%s&%s", url.0, url.1, url.2, url.3, url.4)
   x <- content(GET(url))
   tbls <- lapply(x[[1]], as.data.frame)
   tbl <- do.call(rbind, tbls)
   tbl

} # gtex.demo
#----------------------------------------------------------------------------------------------------


library(RPostgreSQL)
suppressWarnings(
  db.access.test <- try(system("/sbin/ping -c 1 khaleesi", intern=TRUE, ignore.stderr=TRUE)))
if(length(db.access.test) == 0)
   stop("khaleesi database server unavailable")

db <- dbConnect(PostgreSQL(), user="trena", password="trena", dbname="genereg2021", host="khaleesi")
#dbGetQuery(db, "select count(*) from eqtls")
dbGetQuery(db, "select * from eqtls limit 3")
getGenomicRegion(igv)

roi <- getGenomicRegion(igv)
printf("%d", with(roi, end-start))
query <- with(roi, sprintf("select * from eqtls where chrom='%s' and hg19 > %d and hg19 < %d",
                           chrom, start, end))
tbl.eqtl <- dbGetQuery(db, query)
dim(tbl.eqtl)
tbl.eqtl <- subset(tbl.eqtl, pvalue < 0.05)
dim(tbl.eqtl)
tbl.eqtl
tbl.sub <- unique(tbl.eqtl[, c("chrom", "hg19", "rsid", "genesymbol")])

target.genes <- sort(unique(tbl.sub$genesymbol))
for(gene in target.genes){
    tbl.track <- subset(tbl.sub, genesymbol==gene)
    tbl.track$end <- tbl.track$hg19
    tbl.track$start <- tbl.track$hg19 -1
    tbl.track <- tbl.track[, c("chrom", "start", "end", "rsid")]
    track <- DataFrameAnnotationTrack(gene, tbl.track, color="random", trackHeight=25)
    displayTrack(igv, track)
    } # for target.gene

genes.assoc <- c("PPOX","APOA2","ADAMTS4","B4GALT3","USP21","NDUFS2","TOMM40")
intersect(target.genes, genes.assoc)
