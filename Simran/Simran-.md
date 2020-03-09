\#\#Loading Libraries

    library(magrittr)
    library(dplyr)

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    library(tidyverse)

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ stringr 1.4.0
    ## ✓ tidyr   1.0.2     ✓ forcats 0.4.0
    ## ✓ readr   1.3.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x tidyr::extract()   masks magrittr::extract()
    ## x dplyr::filter()    masks stats::filter()
    ## x dplyr::lag()       masks stats::lag()
    ## x purrr::set_names() masks magrittr::set_names()

    library(ggplot2)
    library(reshape2)

    ## 
    ## Attaching package: 'reshape2'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

    library(edgeR)

    ## Warning: package 'edgeR' was built under R version 3.6.1

    ## Loading required package: limma

    library(pheatmap)
    library(ggplot2)

\#Read in Data and Prepare Data Frames

    demo<-read.csv("/Users/simransamra/R/git_temp/Repo_team_Genome-Surfers_W2020/data/raw_data/tcga_paad_clinical.csv", header = T)
    load("/Users/simransamra/R/git_temp/Repo_team_Genome-Surfers_W2020/data/raw_data/tcga_paad.RData")

    #Melt data
    tcga$gene <- rownames(tcga)
    dataMelt<-melt(tcga, id.vars = "gene", var = "Samples")

    #Set up data
    demo <- demo %>% select(c(submitter_id, age_at_index,year_of_birth, year_of_death, vital_status, race, gender, ajcc_pathologic_m, ajcc_pathologic_n,ajcc_pathologic_t, ajcc_pathologic_stage))

    tcgaT <- as.data.frame(t(as.matrix(tcga)))
    tcgaN <- tibble::rownames_to_column(tcgaT, "submitter_id")

    demo$submitter_id <- as.factor(demo$submitter_id)
    data <- right_join(x = tcgaN, y = demo, by = "submitter_id")

    ## Warning: Column `submitter_id` joining character vector and factor, coercing
    ## into character vector

\#Density Plot - Age at Index vs Vital Status

    ggplot(data, aes(x = age_at_index, colour=vital_status)) + 
    geom_density() +
    labs(title = "Age at Index and Vital Status",
    color = "Vital Status", x = "Age at Index", y = "Density") +
    theme(
    plot.title = element_text(color = "blue", size = 12, face = "bold")) +
    scale_fill_manual(values = c("darkblue", "darkred"))

![](Simran-_files/figure-markdown_strict/fig1-1.png) \#Density Plot -
Age at Index vs Pathologic Stage M

    ggplot(data, aes(x = age_at_index, colour=ajcc_pathologic_m)) + 
    geom_density() +
    labs(title = "AJCC Pathologic M",
    subtitle = "AJCC TNM system: (M) Classifies cancers by the presence or absence of distant metastases",  x = "Age at Index", y = "Density", color = "AJCC Pathologic M", caption = "M0:No evidence of distant metastasis; M1: Distant metastasis; MX: Unknown distant metastasis status") +
    theme(
    plot.title = element_text(color = "blue", size = 12, face = "bold"),
    plot.subtitle = element_text(color = "black", size = 7),
    plot.caption = element_text(color = "black", size = 6, hjust = 0)
    )

![](Simran-_files/figure-markdown_strict/fig2-1.png)

\#Density Plot - Age at Index vs Pathologic Stage N

    ggplot(data, aes(x = age_at_index, colour=ajcc_pathologic_n)) + 
    geom_density() +
    labs(title = "AJCC Pathologic N",
    subtitle = "AJCC TNM system: (N) Describes involvement of regional lymph nodes",  x = "Age at Index", y = "Density", color = "AJCC Pathologic N", caption = "--: Not reported; N0: No regional lymph node metastasis; N1: Regional lymph node metastasis; \n  N1b: Metastasis in multiple regional lymph nodes; NX: Metastasis cannot be assessed") +
    theme(
    plot.title = element_text(color = "blue", size = 12, face = "bold"),
    plot.subtitle = element_text(color = "black", size = 7),
    plot.caption = element_text(color = "black", size = 6, hjust = 0)
    )

![](Simran-_files/figure-markdown_strict/fig3-1.png)

\#Density Plot - Age at Index vs Pathologic Stage T

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

![](Simran-_files/figure-markdown_strict/fig4-1.png) ![TMN
Classification](https://www.researchgate.net/publication/279306792/figure/tbl1/AS:601719029379081@1520472404615/TNM-Classification-for-Pancreatic-Cancer-a_W640.jpg)

![Definition of the different pathologic
stages:](https://www.researchgate.net/publication/279306792/figure/tbl2/AS:601719029383177@1520472404649/TNM-Staging-of-Pancreatic-Cancer-a_W640.jpg)

\#Density Plot - Age at Index vs Pathologic Stage

    ggplot(data, aes(x = age_at_index, colour=ajcc_pathologic_stage)) + 
    geom_density() +
    labs(title = "AJCC Pathologic Stage",  x = "Age at Index", y = "Density", color = "AJCC Pathologic Stage") +
    theme(
    plot.title = element_text(color = "blue", size = 12, face = "bold"),
    plot.subtitle = element_text(color = "black", size = 7),
    plot.caption = element_text(color = "black", size = 6, hjust = 0)
    )

![](Simran-_files/figure-markdown_strict/fig5-1.png)

    #Boxplot - Distribution of Gene Expression
    ggplot(dataMelt, aes(x=Samples, y=value)) + 
       geom_boxplot() + 
       xlab("Samples")  + 
      ylab("Expression (Log_2_ Transformed)")+ 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle("Distribution of Gene Expression Across All 177 Samples")+
      theme(plot.title = element_text(hjust = 0.5))

![](Simran-_files/figure-markdown_strict/fig6-1.png)

    #Density Plot
    ggplot(dataMelt, aes(value, color = Samples)) + 
      geom_density() + 
      xlab("Expression (Log2 Transformed)") +
      ylab("Density")+ 
      ggtitle("Distribution of Gene Expression Across All 177 Samples") +
      theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="none")

![](Simran-_files/figure-markdown_strict/fig7-1.png)

    #Prepare Data
    load("/Users/simransamra/R/git_temp/Repo_team_Genome-Surfers_W2020/data/raw_data/tcga_paad.RData")
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
        fontsize = 4)

![](Simran-_files/figure-markdown_strict/fig8-1.png)

    #Prepare Data
    load("/Users/simransamra/R/git_temp/Repo_team_Genome-Surfers_W2020/data/raw_data/tcga_paad.RData")
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
        fontsize = 4)

![](Simran-_files/figure-markdown_strict/fig9-1.png)

    #Prepare Data
    load("/Users/simransamra/R/git_temp/Repo_team_Genome-Surfers_W2020/data/raw_data/tcga_paad.RData")
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
        fontsize = 4)

![](Simran-_files/figure-markdown_strict/fig10-1.png)