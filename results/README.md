# Folder roadmap

The `Oral Presentation` folder contains the .pptx and .docx files as a presentation and a description for it.

Below summarizes the main result where code for a more detailed analysis can be accessed in the `code.Rmd` under the `src` folder.

# Main Results

* The top 5 genes selected using limma, ebayes, and a correlation heatmap is THBS1, NNMT, CREM, ITPRIP, and RP11. The most important variables selected based on LASSO are age_at_index, year_of_diagnosis, race, gender, pathologic T, N, and M, pathologic stage, survival time, THBS1, CREM, and RP11.

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