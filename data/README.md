# **Data**
This directory contains data that was collected from analysis performed with tools from the Andersen Lab, analysis performed with custom computational pipelines utilizing various packages, analysis performed with tools from other publications, and _CaeNDR_ data. 
## **Directory Structure** 
data/ 
- ECs
- HTLDA_analysis
- design_files
- gene_CNV_analysis
- heritability
- filtered_heritability
- target_gene_consequences
## **ECs**
Contains the estimated effective concentrations output files from the custom computational pipelines 
## **HTLDA_analysis**
Contains the raw animal length data extracted from the _ImageXpress_ well images processed using _NemaSize_ for each drug (https://github.com/AndersenLab/NemaSize) 
## **design_files** 
Contains custom design files that layout the details for each drug, metadata plate, and metadata well which are necessary for the _easyXpress_ package (https://github.com/AndersenLab/easyXpress)
## **gene_CNV_analysis** 
Contains
1. gene_names/
   - Contains the identified orthogroups, groups of orthologous genes, of interest
2. orthoFinder
   - Contains the output from _OrthoFinder_ by Emms et al. 2026 (https://github.com/OrthoFinder/OrthoFinder)
## **heritability** 
Contains the broad-sense and narrow-sense heritability calculation output files from the custom computational pipelines 
## **filtered_heritability** 
Contains the broad-sense and narrow-sense heritability calculation output files from the custom computational pipelines using the filtered dataset ( 5 or more wells)
## **proteomes** 
Contains necessary proteomes files for each _Caenorhabditis Briggsae_ strain to run _OrthoFinder_
## **target_gene_consequences** 
Contains the list of variants of implicated genes associated with our strains

