Exploratory Analysis Using All Gene Expression
================

## Loading Libraries

``` r
library(magrittr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("edgeR")
library(edgeR)
library(pheatmap)
library(ggplot2)

# library(devtools)
# install_github("vqv/ggbiplot")
library(ggbiplot)

library(ggpubr)
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(pheatmap))
```

# Read in Data and Prepare Data Frame

``` r
demo<-read.csv("~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/data/raw_data/tcga_paad_clinical.csv", header = T)
load("~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/data/raw_data/tcga_paad.RData")

#Melt data
tcga$gene <- rownames(tcga)
dataMelt<-melt(tcga, id.vars = "gene", var = "Samples")

#Set up data
demo <- demo %>% select(c(submitter_id, age_at_index,year_of_birth, year_of_death, vital_status, race, gender, ajcc_pathologic_m, ajcc_pathologic_n,ajcc_pathologic_t, ajcc_pathologic_stage))

tcgaT <- as.data.frame(t(as.matrix(tcga)))
tcgaN <- tibble::rownames_to_column(tcgaT, "submitter_id")

demo$submitter_id <- as.factor(demo$submitter_id)
data <- right_join(x = tcgaN, y = demo, by = "submitter_id")
```

# Density Plot - Age at Index vs Vital Status

``` r
ggplot(data, aes(x = age_at_index, colour=vital_status)) + 
geom_density() +
labs(title = "Age at Index and Vital Status",
color = "Vital Status", x = "Age at Index", y = "Density") +
theme(
plot.title = element_text(color = "blue", size = 12, face = "bold")) +
scale_fill_manual(values = c("darkblue", "darkred"))
```

![](analysis_SS_files/figure-gfm/fig1-1.png)<!-- -->

# Density Plot - Age at Index vs Pathologic Stage M

``` r
ggplot(data, aes(x = age_at_index, colour=ajcc_pathologic_m)) + 
geom_density() +
labs(title = "AJCC Pathologic M",
subtitle = "AJCC TNM system: (M) Classifies cancers by the presence or absence of distant metastases",  x = "Age at Index", y = "Density", color = "AJCC Pathologic M", caption = "M0:No evidence of distant metastasis; M1: Distant metastasis; MX: Unknown distant metastasis status") +
theme(
plot.title = element_text(color = "blue", size = 12, face = "bold"),
plot.subtitle = element_text(color = "black", size = 7),
plot.caption = element_text(color = "black", size = 6, hjust = 0)
)
```

![](analysis_SS_files/figure-gfm/fig2-1.png)<!-- -->

# Density Plot - Age at Index vs Pathologic Stage N

``` r
ggplot(data, aes(x = age_at_index, colour=ajcc_pathologic_n)) + 
geom_density() +
labs(title = "AJCC Pathologic N",
subtitle = "AJCC TNM system: (N) Describes involvement of regional lymph nodes",  x = "Age at Index", y = "Density", color = "AJCC Pathologic N", caption = "--: Not reported; N0: No regional lymph node metastasis; N1: Regional lymph node metastasis; \n  N1b: Metastasis in multiple regional lymph nodes; NX: Metastasis cannot be assessed") +
theme(
plot.title = element_text(color = "blue", size = 12, face = "bold"),
plot.subtitle = element_text(color = "black", size = 7),
plot.caption = element_text(color = "black", size = 6, hjust = 0)
)
```

![](analysis_SS_files/figure-gfm/fig3-1.png)<!-- -->

# Density Plot - Age at Index vs Pathologic Stage T

``` r
ggplot(data, aes(x = age_at_index, colour=ajcc_pathologic_t)) + 
geom_density() +
labs(title = "AJCC Pathologic T",
subtitle = "AJCC TNM system: (T) ",  x = "Age at Index", y = "Density", color = "AJCC Pathologic T", caption = "--: Not reported; T1: Tumor limited to the pancreas (2 cm or less in greatest dimension); T2: Tumor limited to \n the pancreas (greater than 2 cm in greatest dimension); T3: Tumor extends beyond pancrease, but without \n the involvement of coeliac axis or superior mesenteric artery; T4: Tumor involves coeliac axis or superior \n mesenteric artery; TX: Tumor cannot be assessed )


") +
theme(
plot.title = element_text(color = "blue", size = 12, face = "bold"),
plot.subtitle = element_text(color = "black", size = 7),
plot.caption = element_text(color = "black", size = 6, hjust = 0)
)
```

![](analysis_SS_files/figure-gfm/fig4-1.png)<!-- --> ![TMN
Classification](https://www.researchgate.net/publication/279306792/figure/tbl1/AS:601719029379081@1520472404615/TNM-Classification-for-Pancreatic-Cancer-a_W640.jpg)

![Definition of the different pathologic
stages:](https://www.researchgate.net/publication/279306792/figure/tbl2/AS:601719029383177@1520472404649/TNM-Staging-of-Pancreatic-Cancer-a_W640.jpg)

# Density Plot - Age at Index vs Pathologic Stage

``` r
ggplot(data, aes(x = age_at_index, colour=ajcc_pathologic_stage)) + 
geom_density() +
labs(title = "AJCC Pathologic Stage",  x = "Age at Index", y = "Density", color = "AJCC Pathologic Stage") +
theme(
plot.title = element_text(color = "blue", size = 12, face = "bold"),
plot.subtitle = element_text(color = "black", size = 7),
plot.caption = element_text(color = "black", size = 6, hjust = 0)
)
```

![](analysis_SS_files/figure-gfm/fig5-1.png)<!-- -->

# Heatmaps: Show Correlation

### ajcc pathoogic stages and vital status

``` r
#Prepare Data
load("~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/data/raw_data/tcga_paad.RData")
#tcga$gene <- rownames(tcga)
demoHeat<-demo
toDelete <- seq(1, nrow(demoHeat), 2)
demoHeat<-demoHeat[toDelete ,]

designFactors <- as.data.frame(demoHeat[, c("ajcc_pathologic_stage", "vital_status")])

rownames(designFactors) <- colnames(tcga)
data.matrix <- cor(tcga)
pheatmap(data.matrix, cluster_rows = T, scale = "none", clustering_method = "average", 
    clustering_distance_cols = "correlation", show_colnames = T, show_rownames = T, 
    main = "Clustering Heatmap: Pathologic Stage and Vital Status ", annotation = designFactors, treeheight_col = 35, treeheight_row = 35,
    fontsize = 3)
```

![](analysis_SS_files/figure-gfm/fig8-1.png)<!-- -->

### Race and Gender

``` r
#Prepare Data
load("~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/data/raw_data/tcga_paad.RData")
#tcga$gene <- rownames(tcga)
demoHeat<-demo
toDelete <- seq(1, nrow(demoHeat), 2)
demoHeat<-demoHeat[toDelete ,]

designFactors <- as.data.frame(demoHeat[, c("race",  "gender")])

rownames(designFactors) <- colnames(tcga)
data.matrix <- cor(tcga)
pheatmap(data.matrix, cluster_rows = T, scale = "none", clustering_method = "average", 
    clustering_distance_cols = "correlation", show_colnames = T, show_rownames = T, 
    main = "Clustering Heatmap: Gender and Race", annotation = designFactors, treeheight_col = 35, treeheight_row = 35,
    fontsize = 3)
```

![](analysis_SS_files/figure-gfm/fig9-1.png)<!-- -->

### ajcc pathoogic stages

``` r
#Prepare Data
load("~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/data/raw_data/tcga_paad.RData")
#tcga$gene <- rownames(tcga)
demoHeat<-demo
toDelete <- seq(1, nrow(demoHeat), 2)
demoHeat<-demoHeat[toDelete ,]

designFactors <- as.data.frame(demoHeat[, c("ajcc_pathologic_n", "ajcc_pathologic_t", "ajcc_pathologic_m")])

rownames(designFactors) <- colnames(tcga)
data.matrix <- cor(tcga)
pheatmap(data.matrix, cluster_rows = T, scale = "none", clustering_method = "average", 
    clustering_distance_cols = "correlation", show_colnames = T, show_rownames = T, 
    main = "Clustering Heatmap: AJCC TNM System", annotation = designFactors, treeheight_col = 35, treeheight_row = 35,
    fontsize = 3)
```

![](analysis_SS_files/figure-gfm/fig10-1.png)<!-- -->

# PCA plots for Gener, Race, Vital Status and Pathologic stages

``` r
tcgaT<-t(tcga)
tcga.pca <- prcomp(tcgaT, center = TRUE)
#PCA - Gender
a<-ggbiplot(tcga.pca , var.axes = FALSE, groups=demoHeat$gender)
a + labs(title = "Gender") 
```

![](analysis_SS_files/figure-gfm/fig11%20figure-1.png)<!-- -->

``` r
#PCA - Race
b<-ggbiplot(tcga.pca , var.axes = FALSE, groups=demoHeat$race)
b + labs(title = "Race") 
```

![](analysis_SS_files/figure-gfm/fig11%20figure-2.png)<!-- -->

``` r
#PCA - Vital Status
c<-ggbiplot(tcga.pca , var.axes = FALSE, groups=demoHeat$vital_status)
c + labs(title = "Vital Staus") 
```

![](analysis_SS_files/figure-gfm/fig11%20figure-3.png)<!-- -->

``` r
#PCA - AJCC Pathologic M
d<-ggbiplot(tcga.pca , var.axes = FALSE, groups=demoHeat$ajcc_pathologic_m)
d + labs(title = "AJCC Pathologic M") 
```

![](analysis_SS_files/figure-gfm/fig11%20figure-4.png)<!-- -->

``` r
#PCA - AJCC Pathologic N
e<-ggbiplot(tcga.pca , var.axes = FALSE, groups=demoHeat$ajcc_pathologic_n)
e + labs(title = "AJCC Pathologic N") 
```

![](analysis_SS_files/figure-gfm/fig11%20figure-5.png)<!-- -->

``` r
#PCA - AJCC Pathologic T
f<-ggbiplot(tcga.pca , var.axes = FALSE, groups=demoHeat$ajcc_pathologic_t)
f + labs(title = "AJCC Pathologic T") 
```

![](analysis_SS_files/figure-gfm/fig11%20figure-6.png)<!-- -->

``` r
#PCA - AJCC Pathologic Stage
g<-ggbiplot(tcga.pca , var.axes = FALSE, groups=demoHeat$ajcc_pathologic_stage)
g + labs(title = "AJCC Pathologic Stage") 
```

![](analysis_SS_files/figure-gfm/fig11%20figure-7.png)<!-- -->
