Progress Report
================

# What has changed based on the final proposal (2 pt.)

## Did your dataset change? If so, why?

The datasets did not change. The gene expression dataset was
crossed-referenced with a meta-data containing covariates like age,
gender, survival time, etc. of each subject. Data wrangling was
performed to clean the datasets to facilitate
analysis.

## Have you decided to do a different analysis than what was mentioned in your proposal? If so, Why?

The analysis methods mentioned in the proposal were carried out. Instead
of focusing on the relationship between the gene expression level and
survival time (continuous variable), the goal has been modified to
examine the gene expression level and vital status (dichotomous
variable). The reason is that after a sequence of EDAs, clustering
analysis, and survival analysis, there’s no statistically significant
association between survival time and gene expression of cancer patients
due to the problem of right censoring.

## Are there any changes in task assignments of group members?

No changes in task assignment. All four group members performed a number
of exploratory data analyses. Hassan provided in-context insight and
added a glossary explaining the meaning of variables in the dataset to
help the rest of the group to understand the data better. Hassan and
Simran performed feature selection using PCA, analyzed the correlation
among gene expressions using heatmaps and helped with biological
interpretations. They also helped to organize the repo according to the
rubric’s requirement. Lily and Sophia performed data cleaning, EDA, PCA
and clustering plots, limma with a heat map of the top 10 genes, LASSO,
gene classification, vital status classification, and survival analysis.
Sophia and Lily also wrote the “progress of the analyses” and “results”
sections of the
report.

# What is the progress of the analyses (5 pts.)

## Explain your methodology and progress for the aims you have investigated so far. Which parts were modified and which parts remained the same?

Due to the untidiness of the datasets, data wrangling was an intensely
involved procedure in this analysis. Data cleaning was initially
performed. Subsequently, exploratory data analysis (EDA) was performed
to visualize features and relationships of the covariates like age,
gender, race, and pathological stages using density plots. A heatmap of
correlation among gene expressions was also constructed. However, due to
the presence of a large number of subjects, this heat map is too complex
to infer meaningful conclusions, therefore feature selection on genes
using limma and empirical Bayes was carried out.

To assess differential gene expression, limma and empirical Bayes were
used to determine the top 10 most influential genes on the vital status
of cancer subjects. A heatmap showing correlation among these 10 genes
were constructed, and a further variable selection step was performed.
Within a pair of highly correlated genes, the less significant one
deemed by limma was removed. A cluster dendrogram was plotted to
visualize clustering subjects based on gene expression using the

Euclidean distance between genes. A principal component analysis (PCA)
was conducted to explore the top 10 highly expressed genes from a linear
model fits and to visualize these genes in a lower-dimensional space.
The PCA is able to reduce the gene data dimension by decomposing the
top-ranked genes into different components, and each component is a
distinct subspace that tends to provide the best approximation of the
entire data. The Scree plot shows the proportion of total variance that
was explained by each principal component. To better visualize the
differences in gene expression for patients who were alive compared to
those who did not, we plotted the 10 gene expression values in a 2D
principal component coordinate system in which every gene value is
represented by a new projected (x,y) value.

The least absolute shrinkage and selection operator (LASSO) with 3-fold
cross validation was employed to further reduce the dimension of the
feature space before performing a dichotomous classification (dead or
alive). Variables including age, gender, race, top 5 most important
genes, etc are included as predictors in the LASSO model with vital
status as the outcome. Subsequently, the important features selected
using LASSO were used as the regressors in the random forest and the
Adaptive Boosting (Adaboost) algorithms. The algorithms were trained in
a training set with sample size equaling to the floor value of 75% of
the total sample size. Afterward, the algorithms were run on a
validation set (the remaining 25% of the sample) to obtain accuracy and
area under the receiver operating characteristics curve (known as AUC
for ROC) where AUC ranges from 0 to 1.

Survival analysis was performed using the Cox Proportional Hazard (Cox
PH) model. The covariate in the Cox PH model is `treatment_or_therapy`
containing 3 levels: treatment, therapy, and not reported; it is an
indicator related to the administration of therapeutic agents received.
The goal is to model and to test whether the patients’ survival times
differ among these groups in a statistically significant fashion. A
Kaplan Meier curve was plotted to visualize the survival function. A
log-rank test was performed to assess the significance of a difference
in survival times.

The choice of survival model was modified to Cox PH instead of the
Accelerated Failure Time (AFT) model since the hazard ratio is
approximately constant for all cancer patients and that some patients
immediately died upon diagnosis resulting in survival time (in years) of
zero which further making the AFT model unsuitable. The general roadmap
of the analysis remained the same, namely, data cleaning, EDA,
expression analysis, feature selection, classification, and survival
analysis were carried out as initially
planned.

## What R packages or other tools are you using for your analyses? You do not need to provide your scripts in your report.

R packages used are tidyverse, reshape2, survival, survminer, stringr,
ggpubr, ggridges, glmnet, limma, factoextra, pheatmap, randomForest,
fastAdaboost, and pROC. Tool used is the programming language
R.

## Provide the links to any markdown reports within your repo to refer to the relevant analysis.

[link to code script
(.Rmd)](https://github.com/STAT540-UBC/Repo_team_Genome-Surfers_W2020/blob/master/src/code.Rmd)
[link to knitted code
(.md)](https://github.com/STAT540-UBC/Repo_team_Genome-Surfers_W2020/blob/master/src/code.md)

## Provide references.

Cancer Genome Atlas Research Network. Electronic address:
<andrew_aguirre@dfci.harvard.edu>; Cancer Genome Atlas Research Network.
Integrated Genomic Characterization of Pancreatic Ductal Adenocarcinoma.
Cancer Cell. 2017;32(2):185–203.e13.
<doi:10.1016/j.ccell.2017.07.007>

# Results (3 pts.)

## What are your primary results? Were you able to answer your hypothesis? Did you have any positive results? If no, postulate a discussion as to why that may be. Provide plots and/or tables to present your results. - List some challenges that you encountered? How will you address them?

We were able to answer our hypotheses. Below is a brief summary of the
primary results:

The top 5 genes selected using limma, ebayes, and a correlation heatmap
is THBS1, NNMT, CREM, ITPRIP, and RP11.

![Heatmap of top 10 genes based on expression
level:](~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/docs/heatmap.png)

The most important variables selected based on LASSO are age\_at\_index,
year\_of\_diagnosis, race, gender, pathologic T, N, and M, pathologic
stage, survival time, THBS1, CREM, and RP11.

In the hierarchical dendrogram (clustering) the samples are not
well-separated when we stop the clustering algorithm at k equal to 5. It
might be because of the similarity of gene expressions among all
pancreatic adenocarcinoma patients.

![Hierarchical
dendrogram:](~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/docs/dendro.png)

The Scree plot shows that the first principal component captures 63.1%
of the total variance, and the second principal component explains 8.3%
variance from the remaining variation. We kept the first two principal
components with 71.4% variation retained to compare projected gene
expression points of samples’ vital status. The concentration ellipses
of alive and dead are greatly overlapping, meaning that the gene
expression of alive samples might be similar to that of dead samples.

![Scree
Plot:](~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/docs/scree.png)
![Clustering vital status based on THBS1 and NNMT gene
expression:](~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/docs/clust.png)
![Clustering vital status based on THBS1 and CREM gene
expression:](~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/docs/clust2.png)

![Clustering genes based on principal components 1 and
2:](~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/docs/clust3.png)

![Clustering pathologic stages based on principal components 1 and
2:](~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/docs/clust4.png)

![Clustering vital status based on principal components 1 and
2:](~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/docs/clust5.png)

The accuracy for the random forest algorithm is 70% while that for the
AdaBoost algorithm is 73%. The AUC for the random forest is 0.8825 while
that for AdaBoost is 0.8194. Both algorithms yielded similar results.

![ROC for Random Forest (black) and
AdaBoost(red):](~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/docs/roc.png)
The determining factor of correct classifications is suspected to be
covariates other than gene expression since based on the gene expression
analysis, there is no significant difference in subjects’ top gene
expression. Based on the Kaplan Meier plot, the three group does not
show a significantly different survival rate. A Log-rank test with a
p-value of 0.2 further confirmed that there’s no difference in survival
rate among the groups: treatment, therapy, and not reported.

\[Kaplan Meier Estimate
Curve:\](~/Desktop/git\_docs/Repo\_team\_Genome-Surfers\_W2020/docs/km.png

![Hazard Ratio
Curve:](~/Desktop/git_docs/Repo_team_Genome-Surfers_W2020/docs/hr.png)

One of the primary goals is to identify which target genes affect the
survival rate. However, the real challenge is that we could not find the
target genes from the top 10 genes. One possibility is that the target
genes were not in the top 10 genes in the first place. Secondly, even if
the target genes are located in the group of top 10 genes if a set of
genes are enhancers while the others are inhibitors, the positive and
negative effects that they impose on patients’ survival could be
canceled out. Hence, finding those target genes is an immense challenge
without further data on the genes and patients. The classification
accuracy for both random forest and AdaBoost is semi-optimal. Different
linear combinations of the covariates and different tree depths were
used in both algorithms, but accuracy cannot be improved significantly.
Therefore, an optimization strategy for both algorithms should be
explored. Stratification should be employed to group subjects based on
baseline covariates. However, stratification is not possible with the
current set of covariates because they cannot serve as cluster factors.
More information on subjects is needed. Another challenge is that the
genes may not directly affect survival time, hence making it difficult
to pinpoint target genes solely based on the outcomes of interest -
survival time and vital status.
