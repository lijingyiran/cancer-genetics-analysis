Final Proposal
================

**Motivation and Background**
=============================

Pancreatic adenocarcinoma (PAAD) is the most prevalent form of pancreatic cancer. With a 5 year survival of less than 10%, it is highly aggressive as 80-90% of patients present with surgically unresectable disease \[1-2\]. We hypothesize that certain genes are more heavily expressed in patients with a longer survival time while others seem to associate with shorter survival. Identifying gene targets would enable us to identify ideal chemotherapeutic drug candidates and create a tailored therapy regimen.

**Division of Labor**
=====================

<table style="width:72%;">
<colgroup>
<col width="8%" />
<col width="15%" />
<col width="12%" />
<col width="18%" />
<col width="18%" />
</colgroup>
<thead>
<tr class="header">
<th>Name</th>
<th>Background</th>
<th>Degree</th>
<th>Affiliation</th>
<th>Assigned Task</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Hassan Ali</td>
<td>Biology and CS</td>
<td>Interdisciplinary Oncology</td>
<td>UBC/ BC Cancer / Pancreas Centre BC</td>
<td>Background research and interpretation</td>
</tr>
<tr class="even">
<td>Lily Xia</td>
<td>Stats and CS</td>
<td>Statistics</td>
<td>UBC</td>
<td>Statistical analysis and coding</td>
</tr>
<tr class="odd">
<td>Jingyiran Li</td>
<td>Stats and CS</td>
<td>Statistics</td>
<td>UBC</td>
<td>Statistical analysis and coding</td>
</tr>
<tr class="even">
<td>Simran Samra</td>
<td>Biochemistry and Biology</td>
<td>Experimental Medicine</td>
<td>Heart Lung Innovation / UBC</td>
<td>Background research and interpretation</td>
</tr>
</tbody>
</table>

**Dataset Description**
=======================

This hypothesis will be investigated using the TCGA-PAAD dataset from The Cancer Genome Atlas (TCGA). The dataset contains the expression of 44,084 genes from 177 resected tumors of patients suffering from PAAD. The libraries were sequenced on the Illumina HiSeq 2000 \[3\]. Gene expression levels are measured using RPKM (Reads Per Kilobase of transcript per Million mapped reads) in log10 scale. This dataset additionally contains subject demographics, including, gender, age, race, vital status, year of diagnosis, birth year, death year, and other survival-related information. Since the tumors are resected, it is possible that some patients received neoadjuvant chemotherapy with Gemcitabine. No control subjects are used since this dataset only reflects cancer patient’s genome.

**Aims and Methodology**
========================

The final goal of this project is to identify all genes that are more heavily expressed in PAAD patients that have longer survival time compared to those that have a shorter survival time. This goal will be accomplished using the following methods. Firstly, dimension reduction techniques like LASSO and PCA would be used to reduce the number of covariates. Secondly, clustering algorithms like KNN and EM will be applied to cluster all subjects based on gene expression to identify possible structures that might be present in the data. Then, relevant survival analysis techniques like fitting Cox PH model, AFT model where age will be used to determine whether differences in gene expression levels and survival times are statistically significant. Furthermore, a generalized Linear Mixed-effect (GLM) model will be fitted where random intercept and slope of age for each subject would be taken into account. Finally, training subjects will be classified into a vital status (alive or dead) using appropriate classification techniques, and the classification/prediction model that we get from the training set will be tested on the testing set and validated by k-fold cross-validation. The performance of the classification will be assessed using AUC value for the ROC curve, sensitivity, accuracy, and specificity.

**References**
==============

1.  Schneider G, Siveke JT, Eckel F, Schmid RM. Pancreatic cancer: basic and clinical aspects. Gastroenterology 128(6):1606-1625, 2005.

2.  Conroy T, Bachet JB, Ayav A, Huguet F, Lambert A, Caramella C. Current standards and new innovative approaches for treatment of pancreatic cancer.Eur J Cancer 2016;57:10–22.

3.  Cancer Genome Atlas Research Network. Electronic address: <andrew_aguirre@dfci.harvard.edu>; Cancer Genome Atlas Research Network. Integrated Genomic Characterization of Pancreatic Ductal Adenocarcinoma. Cancer Cell. 2017;32(2):185–203.e13. <doi:10.1016/j.ccell.2017.07.007>
