# Folder roadmap

The `Oral Presentation` folder contains the .pptx and .docx files as a presentation and a description for it.

Below summarizes the main result where code for a more detailed analysis can be accessed in the `code.Rmd` under the `src` folder.

# Main Results

* The top 5 genes selected using limma, ebayes, and a correlation heatmap is `THBS1`, `NNMT`, `CREM`, `ITPRIP`, and `RP11`. The most important variables selected based on LASSO are *age_at_index*, *year_of_diagnosis*, *race*, *gender*, *pathologic T, N, and M*, *pathologic stage*, *survival time*, `THBS1`, `CREM`, and `RP11`.

* We kept the first two principal components with 71.4% variation retained to compare projected gene expression points of samples’ vital status. The concentration ellipses of alive and dead are greatly overlapping, meaning that the gene expression of alive samples might be similar to that of dead samples. 

* For classifying vital status with important covariates and top genes: the accuracy for the random forest algorithm is 70% while that for the AdaBoost algorithm is 73%. The AUC for the random forest is 0.8825 while that for AdaBoost is 0.8194.

* The determining factor of correct classifications is suspected to be covariates other than gene expression since based on the gene expression analysis, there is no significant difference in subjects’ top gene expression. 

* A Log-rank test with a p-value of 0.2 further confirmed that there’s no difference in survival rate among the groups: treatment, therapy, and not reported.

* No survival differences in:
    * No treatment pharm group versus treatment pharm group
    * No treatment pharm versus no treatment radiation
    * Not recorded pharm and radiation groups

* Survival differences in:
    * With treatment pharm and radiation groups
    * No treatment radiation and with treatment radiation groups
    * Not recorded pharm group and with treatment radiation group

* Subtle differences exist between radiation and pharmaceutical treatment. Radiation yields a slightly better survival. Log rank test with a borderline p-value.

# Limitations

* The real challenge is that we could not find the target genes from the top 10 genes, One possibility is that the target genes were not in the top 10 genes in the first place, and positive and negative effects that they impose on patients’ survival could be canceled out. Also, a lack of mutational signature data (`Kras`, `tp53`, `cdkn2a` and `smad4`) further limited our ability to accurately predict target genes. 

* The classification accuracy for both random forest and AdaBoost is semi-optimal. An optimization strategy for both algorithms should be explored.

* Stratification should be employed to group subjects based on baseline covariates. Stratification is not possible with the current set of covariates because they cannot serve as cluster factors.

* Genes may not directly affect survival time, hence making it difficult to pinpoint target genes solely based on the outcomes of interest - survival time and vital status

* In this project, we can only access the data for the cancer patients. Including control cohort might help to bring mutational signature data to the table. 