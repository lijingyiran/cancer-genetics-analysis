code
================

``` r
knitr::opts_chunk$set(echo = TRUE)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(fastAdaboost))
suppressPackageStartupMessages(library(pROC))
```

**Data
Wrangling**

``` r
cli <- read.csv("~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/data/raw_data/tcga_paad_clinical.csv", header = T)
toDelete <- seq(0, nrow(cli), 2)
cli <-  cli[-toDelete, ]
cli1 <- cli %>% select(c(submitter_id, age_at_index, 
                         year_of_birth, year_of_death,year_of_diagnosis, 
                         vital_status, race, gender, ajcc_pathologic_m, 
                         ajcc_pathologic_n, 
                         ajcc_pathologic_t, ajcc_pathologic_stage,
                         treatment_or_therapy, treatment_type))
cli1 <- na.omit(cli1)

load("~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/data/raw_data/tcga_paad.RData")

names(tcga) <- substr(names(tcga), 1, 12)
t <- as.data.frame(t(as.matrix(tcga)))
tcga1 <- tibble::rownames_to_column(t, "submitter_id")


#dat <- right_join(x = tcga1, y = cli1, by = "submitter_id")

cli1$year_of_death <- as.numeric(levels(cli1$year_of_death)[cli1$year_of_death])
```

    ## Warning: NAs introduced by coercion

``` r
cli1$year_of_diagnosis <-as.numeric(levels(cli1$year_of_diagnosis)[cli1$year_of_diagnosis])
```

    ## Warning: NAs introduced by coercion

``` r
cli1$year_of_death[is.na(cli1$year_of_death)] <- 2014
cli1$sur_time <- cli1$year_of_death-cli1$year_of_diagnosis

temp <- data.frame(cbind(as.character(cli$submitter_id), as.character(cli$vital_status)))
names(temp) <- c("submitter_id", "vital_status")
```

# Exploratory Data Analysis

Aim: To visualize features and relationships of the covariates like age,
gender, race, and pathological stages using density plots

## Age distribution across gender

``` r
mu <- cli1 %>% 
  group_by(gender) %>%
  summarise(grp.mean = mean(age_at_index))

ggplot(cli1, aes(x = age_at_index))+ 
  geom_density(aes(fill = gender), alpha = 0.4) +
  geom_vline(aes(xintercept = grp.mean, color = gender),
             data = mu, linetype = "dashed") +
  scale_color_manual(values = c("#868686FF", "#EFC000FF"))+
  scale_fill_manual(values = c("#868686FF", "#EFC000FF"))
```

![](code_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Age distribution across race

## Gene Expression Visualization

``` r
expressionMatrix <- tcga %>% rownames_to_column("gene") %>% as_tibble()
expressionMatrix <- na.omit(expressionMatrix)

expressionMatrix.nogene <- t(scale(t(expressionMatrix[,-1])))

meltedExpressionMatrix <- expressionMatrix %>% melt(id = "gene") 

transformGeneExpressionMatrix <- function(expressionMatrix) {
  expressionMatrix <- expressionMatrix %>%
    as.data.frame() %>% 
    column_to_rownames("gene") %>%
    t() %>% as.data.frame() %>% 
    rownames_to_column("submitter_id") %>% 
    melt(id = "submitter_id") %>% 
    as_tibble() %>% 
    select(submitter_id,
           gene = variable, 
           expression = value)
  return(expressionMatrix)
}

getExpressionForSamples <- function(sampleIds, expressionMatrix) {
  # use gene column as row name
  dataFrame <- expressionMatrix %>% 
    as.data.frame() %>% 
    column_to_rownames("gene")
  return(dataFrame[sampleIds])
}

expressionDataForGene <- transformGeneExpressionMatrix(expressionMatrix)
expressionDataForGene
```

    ## # A tibble: 7,802,868 x 3
    ##    submitter_id gene   expression
    ##    <chr>        <fct>       <dbl>
    ##  1 TCGA-OE-A75W TSPAN6      0.588
    ##  2 TCGA-2J-AABT TSPAN6      0.439
    ##  3 TCGA-IB-7886 TSPAN6      0.875
    ##  4 TCGA-IB-AAUU TSPAN6      0.552
    ##  5 TCGA-2J-AAB6 TSPAN6      0.571
    ##  6 TCGA-LB-A8F3 TSPAN6      0.798
    ##  7 TCGA-HZ-A4BH TSPAN6      0.763
    ##  8 TCGA-IB-7646 TSPAN6      0.917
    ##  9 TCGA-YB-A89D TSPAN6      0.609
    ## 10 TCGA-Z5-AAPL TSPAN6      0.442
    ## # … with 7,802,858 more rows

``` r
expressionDataForGene <- expressionDataForGene %>% left_join(cli1, by = "submitter_id")
```

    ## Warning: Column `submitter_id` joining character vector and factor, coercing
    ## into character vector

``` r
DesMat <- model.matrix(~ vital_status, cli1)
dsFit <- lmFit(expressionMatrix.nogene, DesMat)
ebfit <- eBayes(dsFit)
toptab <- topTable(ebfit)
```

    ## Removing intercept from test coefficients

``` r
important.genes <- expressionMatrix[as.numeric(rownames(toptab)),1]

important.genes.exp <-  tcga1 %>% select(c(submitter_id, important.genes$gene))

meta.with.imp.gene <- right_join(cli1, important.genes.exp, by = "submitter_id")
```

    ## Warning: Column `submitter_id` joining factor and character vector, coercing
    ## into character vector

``` r
#heat map of important genes to show correlation
sampleDists <- as.dist(1-cor(important.genes.exp[,-1]))
sampleDistMatrix <- as.matrix(sampleDists)
dist.rowname <- rownames(sampleDistMatrix)
pheatmap(sampleDistMatrix,cluster_rows = T)
```

![](code_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
#deleted the less important (based on toptable) of 
# the high correlation pairs, ref level is the most imp gene


# CTGF, OLFML2B, CTD-2033D15.2, CYP1B1, KCNE4 deleted.
meta.with.imp.gene.final <- meta.with.imp.gene %>% 
  select(-c("CTGF", "OLFML2B", "CTD-2033D15.2", "CYP1B1", "KCNE4"))
```

**Lasso Selection**

Aim: To further reduce the dimension of the feature space before
performing a dichotomous classification (dead or alive) with least
absolute shrinkage and selection operator (LASSO) with 3-fold cross
validation

``` r
meta.with.imp.gene.final <- na.omit(meta.with.imp.gene.final)
meta.with.imp.gene.final$age_at_index <- as.numeric(meta.with.imp.gene.final$age_at_index)
colnames(meta.with.imp.gene.final)[colnames(meta.with.imp.gene.final) == "RP11-21L23.2"] = "RP11"
smp_size <- floor(0.75 * nrow(meta.with.imp.gene.final))

set.seed(400)
train_ind <- sample(seq_len(nrow(meta.with.imp.gene.final)), size = smp_size)

#The algorithms were trained in a training set with sample size equaling to the floor value of #75% of the total sample size
train <- meta.with.imp.gene.final[train_ind, ]

#The algorithms were run on a validation set (the remaining 25% of the sample) to obtain AUC and #ROC
test <- meta.with.imp.gene.final[-train_ind, ]

yvar <- train$vital_status
temp2 <- train[, - which(names(meta.with.imp.gene.final) %in% c("vital_status", "submitter_id", "year_of_birth", "year_of_death"))]
xvars <- model.matrix(yvar ~ ., data = temp2)

cv.lasso.reg <- cv.glmnet(xvars, yvar, alpha = 1, nfolds = 3, 
                          family = "binomial", measure = "mse", 
                          standardize = T)
best.lam <- cv.lasso.reg$lambda.min
best.lam
```

    ## [1] 0.02809165

``` r
coef(cv.lasso.reg, s = best.lam)
```

    ## 35 x 1 sparse Matrix of class "dgCMatrix"
    ##                                                  1
    ## (Intercept)                           1.382392e+03
    ## (Intercept)                           .           
    ## age_at_index                          1.990851e-02
    ## year_of_diagnosis                    -6.873713e-01
    ## raceblack or african american         .           
    ## racenot reported                      .           
    ## racewhite                            -1.966931e-01
    ## gendermale                           -4.783423e-01
    ## ajcc_pathologic_mM1                   2.346579e-02
    ## ajcc_pathologic_mMX                   .           
    ## ajcc_pathologic_nN0                   .           
    ## ajcc_pathologic_nN1                   .           
    ## ajcc_pathologic_nN1b                  .           
    ## ajcc_pathologic_nNX                   3.650374e-01
    ## ajcc_pathologic_tT1                   .           
    ## ajcc_pathologic_tT2                   .           
    ## ajcc_pathologic_tT3                   1.696760e-01
    ## ajcc_pathologic_tT4                   .           
    ## ajcc_pathologic_tTX                   .           
    ## ajcc_pathologic_stageStage I          .           
    ## ajcc_pathologic_stageStage IA         .           
    ## ajcc_pathologic_stageStage IB         .           
    ## ajcc_pathologic_stageStage IIA        .           
    ## ajcc_pathologic_stageStage IIB        5.006946e-01
    ## ajcc_pathologic_stageStage III        .           
    ## ajcc_pathologic_stageStage IV         2.604733e-03
    ## treatment_or_therapynot reported      7.722387e-01
    ## treatment_or_therapyyes              -3.152969e-01
    ## treatment_typeRadiation Therapy, NOS  .           
    ## sur_time                             -8.018183e-01
    ## THBS1                                 9.536650e-02
    ## NNMT                                  .           
    ## CREM                                 -3.447826e+00
    ## ITPRIP                                .           
    ## RP11                                  3.633714e+00

``` r
plot(cv.lasso.reg)
```

![](code_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
lasso.reg <- glmnet(xvars, factor(yvar), alpha = 1, family = 
                      "binomial", standardize = T)  
plot(lasso.reg, xvar = "lambda", label = T)
```

![](code_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

## Classification of Vital Status

``` r
set.seed(400)
train <- subset(train, select = -c(submitter_id, NNMT, ITPRIP))
test <- subset(test, select = -c(submitter_id,NNMT, ITPRIP))
model1 <- randomForest(vital_status ~ ., data = train, importance = TRUE)
model1
```

    ## 
    ## Call:
    ##  randomForest(formula = vital_status ~ ., data = train, importance = TRUE) 
    ##                Type of random forest: classification
    ##                      Number of trees: 500
    ## No. of variables tried at each split: 4
    ## 
    ##         OOB estimate of  error rate: 31.06%
    ## Confusion matrix:
    ##       Alive Dead class.error
    ## Alive    40   26   0.3939394
    ## Dead     15   51   0.2272727

``` r
# Predicting on train set
predTrain <- predict(model1, train, type = "class")
# Checking classification accuracy
table(predTrain, train$vital_status) 
```

    ##          
    ## predTrain Alive Dead
    ##     Alive    66    0
    ##     Dead      0   66

``` r
predValid <- predict(model1, test, type = "class")
# Checking classification accuracy
mean(predValid == test$vital_status)                    
```

    ## [1] 0.75

``` r
table(predValid,test$vital_status)
```

    ##          
    ## predValid Alive Dead
    ##     Alive    10    3
    ##     Dead      8   23

``` r
# To check important variables
importance(model1)      
```

    ##                            Alive       Dead MeanDecreaseAccuracy
    ## age_at_index           3.7837874 -0.3650419           2.60286726
    ## year_of_birth          2.0052041 -2.9793274          -0.26802418
    ## year_of_death         25.2372559 19.1604012          25.19877841
    ## year_of_diagnosis      3.1203853  1.7614282           3.61398131
    ## race                  -1.7120786 -2.4464529          -2.77191367
    ## gender                 3.1957658  5.1286842           5.25012841
    ## ajcc_pathologic_m      0.2257410 -0.8869030          -0.52637393
    ## ajcc_pathologic_n      2.3724764  0.5810740           2.25894283
    ## ajcc_pathologic_t      1.4243925  3.4677205           3.31915383
    ## ajcc_pathologic_stage  3.5063169  3.7966839           5.07213714
    ## treatment_or_therapy   1.7546826  0.9586922           1.78366842
    ## treatment_type        -1.4538082  0.4516831          -1.03927503
    ## sur_time               9.8296904  6.6201143          11.42102577
    ## THBS1                  1.0959312  2.8368672           2.64924974
    ## CREM                  -0.4516851  1.0359117          -0.01997821
    ## RP11                   6.7959284  2.7152865           6.56260315
    ##                       MeanDecreaseGini
    ## age_at_index                 2.6792137
    ## year_of_birth               18.7066015
    ## year_of_death               12.7840591
    ## year_of_diagnosis            2.0546732
    ## race                         0.9167157
    ## gender                       1.4924602
    ## ajcc_pathologic_m            1.0017132
    ## ajcc_pathologic_n            1.0684579
    ## ajcc_pathologic_t            1.1491863
    ## ajcc_pathologic_stage        2.6155060
    ## treatment_or_therapy         1.6359480
    ## treatment_type               0.5911329
    ## sur_time                     5.5470424
    ## THBS1                        3.9562686
    ## CREM                         3.7170173
    ## RP11                         5.4657492

``` r
varImpPlot(model1)
```

![](code_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
#Adaboost
ada <- adaboost(vital_status ~., data =  train, 10, tree_depth = 5, n_rounds = 200, verbose = TRUE)
pred <- predict(ada,newdata=test)
print(pred$error)
```

    ## [1] 0.2954545

``` r
print(table(pred$class,test$vital_status))
```

    ##        
    ##         Alive Dead
    ##   Alive    13    8
    ##   Dead      5   18

``` r
predValid1 <- predict(model1, test, type = "prob")
rf.roc<-roc(test$vital_status,predValid1[,2])
```

    ## Setting levels: control = Alive, case = Dead

    ## Setting direction: controls < cases

``` r
ada.roc<-roc(test$vital_status,pred$votes[,2])
```

    ## Setting levels: control = Alive, case = Dead
    ## Setting direction: controls < cases

``` r
plot(rf.roc)
lines(ada.roc, col = "red")
```

![](code_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
auc(rf.roc)
```

    ## Area under the curve: 0.8846

``` r
auc(ada.roc)
```

    ## Area under the curve: 0.8077

**Classification of Gene Expressions**

## Data set

``` r
# top 10 genes
selectGene <- tcga1 %>% 
  select(c(submitter_id, important.genes$gene))

temp1 <- as.data.frame(selectGene %>% 
                         melt(id = "submitter_id",
                              var = "gene"))

temp2 <- left_join(temp1, cli1) 
```

    ## Joining, by = "submitter_id"

    ## Warning: Column `submitter_id` joining character vector and factor, coercing
    ## into character vector

``` r
first10 <- temp2 %>% 
  pivot_wider(names_from = "gene",
              values_from = "value")

disdat <- first10[,16:25]
# top 10 Dendogram
dist.euclidean <- as.dist(1-cor(disdat))
p <- hclust(dist.euclidean, method = "single")
plot(p, xlim = 50)
rect.hclust(p, k =5)
```

![](code_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
first10 %>% 
  select("submitter_id", "CREM", "THBS1", "vital_status") %>% 
  ggplot(aes(x = CREM, y = THBS1, color = vital_status)) +
  geom_point() +
  theme_bw()
```

![](code_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
tempp <- left_join(temp1, cli) 
```

    ## Joining, by = "submitter_id"

    ## Warning: Column `submitter_id` joining character vector and factor, coercing
    ## into character vector

``` r
first <- tempp %>% 
  pivot_wider(names_from = "gene",
              values_from = "value")
```

## PCA

``` r
res.pca1 <- prcomp(first[,162:171], scale = T)
fviz_eig(res.pca1,addlabels = TRUE)
```

![](code_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
var <- get_pca_var(res.pca1)

set.seed(123)
res.km <- kmeans(var$coord, centers = 2, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(res.pca1, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF"),
             legend.title = "Cluster",
             title = "")
```

![](code_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
fviz_pca_ind(res.pca1,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = factor(first10$ajcc_pathologic_stage), # color by groups
             palette = c("#0073C2FF", "#EFC000FF", "green", "blue",
                         "lightblue", "red", "grey", "black"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "8 ajcc pathologic stages",
             alpha.ind = 0.7,
             title = ""
)
```

![](code_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

``` r
fviz_pca_ind(res.pca1,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = factor(first10$vital_status), # color by groups
             palette = c("#0073C2FF","green"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "vital status",
             alpha.ind = 0.7,
             title = ""
)
```

![](code_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

## Comparing alive and dead based on average top10 gene expresion

``` r
first$count <- scale(rowMeans(first[,162:171]))


first %>% 
  ggplot(aes(x = vital_status, y = count)) +
  geom_violin(width = 0.3) +
  theme_bw()
```

![](code_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
# two sample t-test
t.test(first[first$vital_status == "Alive", ]$count, first[first$vital_status == "Dead", ]$count)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  first[first$vital_status == "Alive", ]$count and first[first$vital_status == "Dead", ]$count
    ## t = -0.29897, df = 164.9, p-value = 0.7653
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.3452312  0.2544316
    ## sample estimates:
    ##   mean of x   mean of y 
    ## -0.02359762  0.02180215

``` r
sample_to_choose <- sample(1:length(unique(first10$submitter_id)), size = 100)
names_to_choose <- as.character(unique(first10$submitter_id)[sample_to_choose])

temp10 <- first %>% 
  filter(submitter_id %in% names_to_choose) %>% 
  group_by(submitter_id) 

temp10 %>% 
  ggplot(aes(x = as.factor(submitter_id), y = count, color = vital_status)) + 
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

![](code_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

**Survival Analysis**

## Kalpan-Meier estimator

We first use Kalpan-Meier plot to summarize the survival experience of
the event-time porcess. `sur_time` variable is the duration from year of
diagnosis to
death.

``` r
sur_dat <- as.data.frame(cbind(cli1$sur_time, cli1$vital_status, cli1$treatment_or_therapy, cli1$treatment_type))
colnames(sur_dat) <- c("sur_time", "vital_status", "treatment_or_therapy", "treatment_type")
with(cli1, Surv(sur_time, vital_status))
```

    ##   [1]  0:Dead  3+      0:Dead  5+      0:Dead  1+      1+      1+      0:Dead
    ##  [10]  1+      2+      6+      1+      1+      1:Dead  2+      0:Dead  1+    
    ##  [19]  2:Dead  3+      1+      1:Dead  3:Dead  1:Dead  4+      1:Dead  1:Dead
    ##  [28]  1+      3:Dead  2+      4+      2:Dead  1+      1+      2+      1:Dead
    ##  [37]  2:Dead  0:Dead  0:Dead  4:Dead  0:Dead  0:Dead  1+      3+      1+    
    ##  [46]  1:Dead  2+      0:Dead  5+      0:Dead  0:Dead  3+      0:Dead  0:Dead
    ##  [55]  1+      3:Dead  4+      2+      6+      1:Dead  4+      2+      0:Dead
    ##  [64]  2:Dead  3+      2+      1:Dead  6+      2+      1:Dead  3+      2:Dead
    ##  [73]  2:Dead  2+      1:Dead  1:Dead  3+      1:Dead  1:Dead  1+      2+    
    ##  [82]  3+      3+      1:Dead  0:Dead  1:Dead  1:Dead  1+      1:Dead  1+    
    ##  [91]  0:Dead  3+      2:Dead  2+      5+      2+      0:Dead  1:Dead  1:Dead
    ## [100]  2:Dead  1:Dead  0:Dead  3:Dead  1+      1:Dead  1:Dead  1:Dead  4:Dead
    ## [109]  3+      0:Dead  0:Dead  3+     NA+     13:Dead  5+      5+      1+    
    ## [118]  1+      3+      1+      3:Dead  2+      2+      3:Dead  1:Dead  0:Dead
    ## [127]  2+      1:Dead  2:Dead  2+      2:Dead  0:Dead  2+      1:Dead  1+    
    ## [136]  1+      1:Dead  7+      0:Dead  2+      3+      2:Dead  3:Dead  2:Dead
    ## [145]  2:Dead  1:Dead  4+      2+      3+      3:Dead  2:Dead  1:Dead  1+    
    ## [154]  1:Dead  1+      5:Dead  1+      2+      2+      2+      0:Dead  0:Dead
    ## [163]  2+      1:Dead  3:Dead  2:Dead  3+      0:Dead  2:Dead  1+      6+    
    ## [172]  1:Dead  3:Dead  2+      1+      3:Dead  0:Dead

``` r
fit1 <- survfit(Surv(sur_time, vital_status) ~ treatment_or_therapy + treatment_type, data = sur_dat)
summary(fit1)
```

    ## Call: survfit(formula = Surv(sur_time, vital_status) ~ treatment_or_therapy + 
    ##     treatment_type, data = sur_dat)
    ## 
    ## 1 observation deleted due to missingness 
    ##                 treatment_or_therapy=1, treatment_type=1 
    ##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
    ##     0     14       3    0.786   0.110        0.598        1.000
    ##     1     11       3    0.571   0.132        0.363        0.899
    ##     2      7       1    0.490   0.136        0.284        0.845
    ## 
    ##                 treatment_or_therapy=1, treatment_type=2 
    ##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
    ##     0     62      18    0.710  0.0576        0.605        0.832
    ##     1     44      12    0.516  0.0635        0.406        0.657
    ##     2     23       4    0.426  0.0664        0.314        0.579
    ##     3     13       3    0.328  0.0714        0.214        0.502
    ##    13      1       1    0.000     NaN           NA           NA
    ## 
    ##                 treatment_or_therapy=2, treatment_type=1 
    ##         time       n.risk      n.event     survival      std.err lower 95% CI 
    ##        1.000        5.000        2.000        0.600        0.219        0.293 
    ## upper 95% CI 
    ##        1.000 
    ## 
    ##                 treatment_or_therapy=2, treatment_type=2 
    ##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
    ##     2      5       1      0.8   0.179        0.516            1
    ##     3      4       2      0.4   0.219        0.137            1
    ##     5      1       1      0.0     NaN           NA           NA
    ## 
    ##                 treatment_or_therapy=3, treatment_type=1 
    ##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
    ##     0     69       7    0.899  0.0363        0.830        0.973
    ##     1     62      14    0.696  0.0554        0.595        0.813
    ##     2     34       7    0.552  0.0653        0.438        0.696
    ##     3     15       5    0.368  0.0801        0.240        0.564
    ##     4      5       1    0.295  0.0919        0.160        0.543
    ## 
    ##                 treatment_or_therapy=3, treatment_type=2 
    ##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
    ##     1     20       2    0.900  0.0671        0.778        1.000
    ##     2     17       3    0.741  0.0999        0.569        0.965
    ##     3      7       1    0.635  0.1302        0.425        0.949
    ##     4      3       1    0.424  0.1935        0.173        1.000

``` r
ggsurvplot(
  fit1, 
  data =sur_dat , 
  size = 1,                 # change line size
  palette = 
    c("green", "dark green", "light blue", "blue","pink","red"),# custom color palettes
  conf.int = F,          # Add confidence interval
  ggtheme = theme_bw(),
  legend = "right",
  legend.title = "treatment levels",
  legend.labs = c("no treatment (Pharmaceutical group)",
                  "no treatment (Radiation group)",
                  "not recorded (Pharmaceutical group)",
                  "not recorded (Radiation group)",
                  "with treatment (Pharmaceutical group)",
                  "with treatment (Radiation group)")
) 
```

![](code_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
# with treatment 
sur_dat2 <- sur_dat %>% filter(treatment_type == 2 & treatment_or_therapy %in% c(1,3)) 
colnames(sur_dat2) <- c("sur_time", "vital_status", "treatment_or_therapy", "treatment_type")
fit2 <- survfit(Surv(sur_time, vital_status) ~ treatment_or_therapy + treatment_type, data = sur_dat2)
summary(fit2)
```

    ## Call: survfit(formula = Surv(sur_time, vital_status) ~ treatment_or_therapy + 
    ##     treatment_type, data = sur_dat2)
    ## 
    ## 1 observation deleted due to missingness 
    ##                 treatment_or_therapy=1, treatment_type=2 
    ##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
    ##     0     62      18    0.710  0.0576        0.605        0.832
    ##     1     44      12    0.516  0.0635        0.406        0.657
    ##     2     23       4    0.426  0.0664        0.314        0.579
    ##     3     13       3    0.328  0.0714        0.214        0.502
    ##    13      1       1    0.000     NaN           NA           NA
    ## 
    ##                 treatment_or_therapy=3, treatment_type=2 
    ##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
    ##     1     20       2    0.900  0.0671        0.778        1.000
    ##     2     17       3    0.741  0.0999        0.569        0.965
    ##     3      7       1    0.635  0.1302        0.425        0.949
    ##     4      3       1    0.424  0.1935        0.173        1.000

``` r
ggsurvplot(
  fit2, 
  data =sur_dat2 , 
  size = 1,                 # change line size
  palette = 
    c("light blue", "blue"),# custom color palettes
  conf.int = T,          # Add confidence interval
  pval = TRUE,              # Add p-value
  ggtheme = theme_bw(),
  legend = "top",
  legend.title = "Radiation treatment",
  legend.labs = c("no treatment",
                  "with treatment")
)
```

![](code_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
# with treatment compare treatment types
sur_dat3 <- sur_dat %>% filter(treatment_or_therapy == 3)
colnames(sur_dat3) <- c("sur_time", "vital_status", "treatment_or_therapy", "treatment_type")
fit3 <- survfit(Surv(sur_time, vital_status) ~ treatment_or_therapy + treatment_type, data = sur_dat3)
summary(fit3)
```

    ## Call: survfit(formula = Surv(sur_time, vital_status) ~ treatment_or_therapy + 
    ##     treatment_type, data = sur_dat3)
    ## 
    ## 1 observation deleted due to missingness 
    ##                 treatment_or_therapy=3, treatment_type=1 
    ##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
    ##     0     69       7    0.899  0.0363        0.830        0.973
    ##     1     62      14    0.696  0.0554        0.595        0.813
    ##     2     34       7    0.552  0.0653        0.438        0.696
    ##     3     15       5    0.368  0.0801        0.240        0.564
    ##     4      5       1    0.295  0.0919        0.160        0.543
    ## 
    ##                 treatment_or_therapy=3, treatment_type=2 
    ##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
    ##     1     20       2    0.900  0.0671        0.778        1.000
    ##     2     17       3    0.741  0.0999        0.569        0.965
    ##     3      7       1    0.635  0.1302        0.425        0.949
    ##     4      3       1    0.424  0.1935        0.173        1.000

``` r
ggsurvplot(
  fit3, 
  data =sur_dat3 , 
  size = 1,                 # change line size
  palette = c("orange","red"),# custom color palettes
  conf.int = T,          # Add confidence interval
  pval = TRUE,              # Add p-value
  ggtheme = theme_bw(),
  legend = "top",
  legend.title = "with treatment",
  legend.labs = c("Pharmaceutical treatment",
                  "Radiation treatment")
)
```

![](code_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

``` r
fit.coxph <- coxph(Surv(sur_time, vital_status) ~ treatment_or_therapy, 
                   data = sur_dat)

ggforest(fit.coxph, data = sur_dat)
```

![](code_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

## Log rank test

``` r
survdiff(Surv(sur_time, vital_status) ~ treatment_or_therapy, 
         data = sur_dat)
```

    ## Call:
    ## survdiff(formula = Surv(sur_time, vital_status) ~ treatment_or_therapy, 
    ##     data = sur_dat)
    ## 
    ## n=176, 1 observation deleted due to missingness.
    ## 
    ##                         N Observed Expected (O-E)^2/E (O-E)^2/V
    ## treatment_or_therapy=1 76       45    36.97    1.7452    3.6620
    ## treatment_or_therapy=2 11        6     6.69    0.0712    0.0951
    ## treatment_or_therapy=3 89       41    48.34    1.1150    2.9390
    ## 
    ##  Chisq= 3.7  on 2 degrees of freedom, p= 0.2

``` r
#p=0.2 no difference in survival
```
