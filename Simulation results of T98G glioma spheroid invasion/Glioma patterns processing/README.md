## Glioma patterns processing ##
This folder contains the code which derives the numerical datasets from the glioma invasion patterns (in binary form). 

1. Each observed glioma pattern (**pattern glioma t1 -40um 400_400 binary.tif** and **pattern glioma t3 -40um 400_400 binary.tif**) was dissected into 100 sections (10x10) using
ImageJ software, the dissected sections are stored in folder **Glioma pattern d1** and **Glioma pattern d3**. 
2. The density of each section was measured by the ratio of black areas in it, the numerical values are stored in the .csv files in the folder. 
3. At the end, file **Glioma data processing.R** reads in the numerical values stored in these .csv files and write them into a single .rds file **Glioma simulation full reference densities.rds**. 
