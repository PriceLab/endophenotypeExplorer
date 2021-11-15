#----------------------------------------------------------------------------------------------------
get.rna.matrix.codes <- function()
{
  list("max-tcx"="mayo.tcx.robinson.normalized.PMI-age-cellType-covariate-collected.15201x262.RData",
       "sage-eqtl-cer"="mtx.mayo.cer.eqtl-optimized-geneSymbols-sampleIDs-17009x261.RData",
       "sage-eqtl-tcx"="mtx.mayo.tcx.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x257.RData",
       "old-mayo-tcx"="mayo.tcx.16969x262.covariateCorrection.log+scale.RData",
       "old-mayo-cer"="cerebellum.15167x263.RData",
       "old-rosmap"="rosmap.14235x632.RData",
       "sage-eqtl-rosmap"="mtx.rosmap.rnaseq-residual-eqtl-geneSymbols-patients-15582x632.RData",
       "sage-counts-rosmap"="mtx.rosmap.rnaseq-counts-geneSymbols-patients-15582x632.RData",
       "max-rosmap"="ROSMAP_rnaseq-ncsNormalized-covCorrectedForStudyPMIandAgeAtDeath-10oct21-21667x631.RData",
       "tcx-unknown"="temporalCortex.15167x264.RData")


} # get.rna.matrix.codes
#----------------------------------------------------------------------------------------------------
get.rna.matrix <- function(code)
{
   mtx.choices <- c("max-tcx", "sage-eqtl-cer", "sage-eqtl-tcx",
                    "old-mayo-tcx", "old-mayo-cer",
                    "old-rosmap", "sage-eqtl-rosmap",
                    "sage-counts-rosmap", "max-rosmap")

   stopifnot(Sys.info()[["user"]] %in% c("paul", "pshannon"))
   stopifnot(code %in% names(get.rna.matrix.codes()))

   if(code=="max-rosmap"){
     data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression/maxNormalizations"
     filename <- "ROSMAP_rnaseq-ncsNormalized-covCorrectedForStudyPMIandAgeAtDeath-10oct21-21667x631.RData"
     mtx.rna <- get(load(file.path(data.dir, filename)))
     }

   if(code=="max-tcx"){
     data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression/maxNormalizations"
     filename <- "mayo.tcx.robinson.normalized.PMI-age-cellType-covariate-collected.15201x262.RData"
     mtx.rna <- get(load(file.path(data.dir, filename)))
     }

   if(code=="sage-eqtl-cer"){
     data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression/sage.eqtl.optimized"
     filename <- "mtx.mayo.cer.eqtl-optimized-geneSymbols-sampleIDs-17009x261.RData"
     mtx.rna <- get(load(file.path(data.dir, filename)))
     }

   if(code=="sage-eqtl-tcx"){
     data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression/sage.eqtl.optimized"
     filename <- "mtx.mayo.tcx.eqtl-optimized-geneSymbols-sampleIDs-with-vcf17009x257.RData"
     mtx.rna <- get(load(file.path(data.dir, filename)))
     }

   if(code=="old-mayo-tcx"){
     data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression"
     filename <- "mayo.tcx.16969x262.covariateCorrection.log+scale.RData"
     mtx.rna <- get(load(file.path(data.dir, filename)))
     }

   if(code=="old-mayo-cer"){
     data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression"
     filename <- "cerebellum.15167x263.RData"
     mtx.rna <- get(load(file.path(data.dir, filename)))
     }

   if(code=="old-rosmap"){
     data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression"
     filename <- "rosmap.14235x632.RData"
     mtx.rna <- get(load(file.path(data.dir, filename)))
     }

   if(code=="sage-eqtl-rosmap"){
     data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression/sage.eqtl.optimized"
     filename <- "mtx.rosmap.rnaseq-residual-eqtl-geneSymbols-patients-15582x632.RData"
     mtx.rna <- get(load(file.path(data.dir, filename)))
     }
   if(code=="sage-counts-rosmap"){
     data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression/sage.eqtl.optimized"
     filename <- "mtx.rosmap.rnaseq-counts-geneSymbols-patients-15582x632.RData"
     mtx.rna <- get(load(file.path(data.dir, filename)))
     }

   if(code=="tcx-unknown"){
     data.dir <- "~/github/TrenaProjectAD/inst/extdata/expression"
     filename <- "temporalCortex.15167x264.RData"
     mtx.rna <- get(load(file.path(data.dir, filename)))
     }

   invisible(mtx.rna)

} # get.rna.matrix
#----------------------------------------------------------------------------------------------------
