This folder comntains all the data files used in this project

### **Raw Data**

The dataset contains the expression of 44,084 genes from 177 resected tumors of patients suffering from PAAD. The libraries were sequenced on the Illumina HiSeq 2000. Gene expression levels are measured using RPKM (Reads Per Kilobase of transcript per Million mapped reads) in log10 scale. This dataset additionally contains subject demographics, including, gender, age, race, vital status, year of diagnosis, birth year, death year, and other survival-related information. Since the tumors are resected, it is possible that some patients received neoadjuvant chemotherapy with Gemcitabine or radiation. No control subjects are used since this dataset only reflects cancer patient’s genome.


### **Processed Data**

* Sample clinical data had multiple unfilled columns, and they were removed. 
* We mainly focused on patients’ age, sex, race, tumor stages, treatments, survival time and vital status.
* To allow for comparability of the effects of these on gene expression, expression data was standardized to generate heatmaps


