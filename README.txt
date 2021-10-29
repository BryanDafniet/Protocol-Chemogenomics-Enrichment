Two directories for the creation of the chemogenomics library and the enrichment method used on the compounds in the database.

For the Chemogenomics_Protocol :

The R script is where the whole protocol is written as well as the NEO4j commands to generate intermediate files out of our database.

The FinalData_5k100 is the file containing the final library.

All the other csv files are the intermediate files generated from R or NEO4j required to run the script and compare our result from an eventual reproduction of our method.

For the Enrichment_Protocol :

The The R script is where the whole protocol is written.

The csv scripts are the various enrichment matrices calculated from the R packages included in our database. 


IMPORTANT NOTE:

2 files needed in the Chemogenomics_Protocol are on a drive due to their size :
	-Dist_Mat_Mol41k_Redone.csv 
	-Res_5000_1000_600.csv

They are accessible here:

https://drive.google.com/drive/folders/1Y0niK56IeFt2Y4U0UWqE37dme0BnybLW?usp=sharing