library(motifbreakR)
print(load("ptk2b.eqtls.top.14.07Sep.0640.motifbreaks.RData"))

x <- system.time(
results.pval <- calculatePvalue(results)
)
print(x)
save(results.pval, file="results.top14.pval.RData")
