STAT 540 Statistical Methods in High Dimensional Biology
=================================================

An Analysis on Gene Expressions of Pancreatic adenocarcinoma (PAAD)
--------------------------------------------------------------------------------
**Repo Navigation:**

- *[Group Information](#org0)*
- *[Project Backgroud](#org1)*
- *[Project Objectives](#org2)*
- *[Directory Roadmap](#org3)*
- *[Main Results](#org4)*



<a id="org0"></a>
Group Name: Genome Surfers

| Team Memembers  |  Reviewers |
| :---: |  :---: |  
| Jingyiran Li |  Sara Mostafavi |
| Lily Xia |  Sina Jafarzadeh |
|  Hassan  |  Victor Yuan |
| Simran Samra |  Keegan Korthauer |

<a id="org1"></a>
### Project background

Pancreatic adenocarcinoma (PAAD) is the most prevalent form of pancreatic cancer. With a 5 year survival of less than 10%, it is highly aggressive as 80-90% of patients present with surgically unresectable disease [1-2]. We hypothesize that certain genes are more heavily expressed in patients with a longer survival time while others seem to associate with shorter survival. Identifying gene targets would enable us to identify ideal chemotherapeutic drug candidates and create a tailored therapy regimen. 


<a id="org2"></a>
### Project objectives

-   To identify all genes that are more heavily expressed in PAAD patients that have longer survival time compared to those that have a shorter survival time.
-   Training subjects will be classified into a vital status (alive or dead) using appropriate classification techniques


<a id="org3"></a>
### Directory roadmap

(We may want to edit this once we have better-defined folders. Below an example.) Each directory includes its own README file. In general: 
* [`docs`](https://github.com/STAT540-UBC/Repo_team_Genome-Surfers_W2020/tree/master/docs) includes the written portions of the project: main, sections, appendices, and project outline tex/pdfs/aux files. 
  + [:point_right: project final proposal](https://github.com/STAT540-UBC/Repo_team_Genome-Surfers_W2020/blob/master/docs/project_proposal.md)
  + [:point_right: progress report](https://github.com/STAT540-UBC/Repo_team_Genome-Surfers_W2020/blob/master/docs/progress_report.md)
  + [:point_right: oral presentation files](https://github.com/STAT540-UBC/Repo_team_Genome-Surfers_W2020/tree/master/results/Oral%20presentation)
* [`src`](https://github.com/STAT540-UBC/Repo_team_Genome-Surfers_W2020/tree/master/src) includes the R code developed for the project.
* [`data`](https://github.com/STAT540-UBC/Repo_team_Genome-Surfers_W2020/tree/master/data) contains data files used for the project.


### Main Results

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


