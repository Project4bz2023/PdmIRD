# PdmIRD version 1.0 #
## About

PdmIRD is a tool to predict IRD (Inherited retinal diseases) mutations. It is applicable to missense mutations in IRD-related genes. We also provide a website (https://zbshiny.shinyapps.io/IRDmis/) for online queries, which integrates the information from the Genome Aggregation Database (gnomAD), the Exome Aggregation Consortium (ExAC), the 1000 genomes project (1KGP), the China Metabolic Analytics Project (ChinaMAP) and the Westlake BioBank for Chinese (WBBC).

## Dependence
### Tested environment:
Perl (version 5.16.3)  
R (version 4.1.0)  
### Required packages:
Boruta_8.0.0         
randomForest_4.7-1.1  
missForest_1.5        
caret_6.0-94          
ggfortify_0.4.16   
pROC_1.18.4   
ggplot2_3.4.3   
tidyverse_2.0.0   
Hmisc_5.1-1         
## How to use
First, using the ANNOVAR (https://annovar.openbioinformatics.org/en/latest/user-guide/download/) to annotate the mutations, the 'step1annovar.sh' in the script folder is an example.
Then, 'step2run.sh' uses the 'annovar.matchOtheranno.pl' to annotate other information. Part of the source data is stored in the data folder, others you can download from the office sites.
Finally, you can obtain the prediction model by executing the 'preprocess_train_test.R', and use it to get your prediction result. 

Of course, you also can use the pre-calculated score file, "IRDmis.predict.gz".
