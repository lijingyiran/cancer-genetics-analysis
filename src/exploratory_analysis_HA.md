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
tcga_qcd <- tcga[, !(colnames(tcga) %in% outlier_samples)]

#PCA of filtered dataset
A <- cov(tcga_qcd, method = "pearson")
E <- eigen(A)
FV <- data.frame(as.matrix(E$vectors[, 1:2]), stringsAsFactors = F)
colnames(FV) <- c("PC1", "PC2")
rownames(FV) <- colnames(tcga_z_qcd)

#Make Sample ID's Comparable
outlier_samples <- substring(outlier_samples, 1, 12)
clean_meta_qcd <- clean_meta[!(clean_meta$submitter_id %in% outlier_samples), ]
rownames(FV) <- substring(rownames(FV),1,12)
```

*Plots*

``` r
#grouped by race
ggplot(FV, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color=clean_meta_qcd$race))
```

![](exploratory_analysis_HA_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#grouped by vital status
ggplot(FV, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color=clean_meta_qcd$vital_status))
```

![](exploratory_analysis_HA_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
#grouped by gender
ggplot(FV, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color=clean_meta_qcd$gender))
```

![](exploratory_analysis_HA_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
#grouped by pathologic stage
ggplot(FV, aes(x=PC1, y=PC2)) + 
  geom_point(aes(color=clean_meta_qcd$vital_status)) +
  facet_wrap(~clean_meta_qcd$ajcc_pathologic_stage)
```

![](exploratory_analysis_HA_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->
