exploratory\_analysis
================

*PCA of top 5000 most variable genes*

``` r
#PCA for top 5000 most variable genes

tcga_z <- scale(tcga)
tcga_var <- as.matrix(as.data.frame(apply(tcga, 1, var)))
tcga_var_sorted <- tcga_var[order(tcga_var, decreasing = TRUE), , drop = FALSE] 
tcga_var_top5000 <- as.matrix(as.data.frame(tcga_var_sorted[1:5000,]))
tcga_filt <- tcga_z[rownames(tcga_z) %in% rownames(tcga_var_top5000),]
A <- cov(tcga_filt, method = "pearson")
E <- eigen(A)
FV <- data.frame(as.matrix(E$vectors[, 1:2]), stringsAsFactors = F)
colnames(FV) <- c("PC1", "PC2")
rownames(FV) <- colnames(tcga_filt)
plot(FV[, 1:2])
```

![](exploratory_analysis_HA_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->
*Cleaning the Dataset*

``` r
#Order by PC1
FV_sorted <- FV[order(-FV$PC1),]

#identify 8 outlier samples
outlier_samples <- rownames(FV_sorted[1:8,])

#Remove 8 outlier samples
tcga_z_qcd <- tcga_z[, !(colnames(tcga_z) %in% outlier_samples)]
```
