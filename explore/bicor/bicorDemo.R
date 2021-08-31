library(WGCNA)

x <- mtcars$disp
y <- mtcars$mpg
WGCNA::bicor(x, y, use="pairwise.complete.obs")

cor(x, y, use="pairwise.complete.obs", method="spearman")
cor(x, y, use="pairwise.complete.obs", method="pearson")
cor(x, y, use="pairwise.complete.obs", method="kendall")
