# **Scripts**
This directory contains the scripts used to carry out analysis and generate figures and tables for the manuscript.
- `EC_calculations.R`
Takes the processed and cleaned animal length data and calculates the estimated effective concentrations (10, 50, and 90), performs effective concentration comparisons between strains within each drug, 
calculates the slope of Dose Response Curves, performs slope comparisons between strains within each drug, calculates the heritability, and calculates the heritability at each effective concentration. 
- `cnv_analysis.R`
Performs the Call Number Variant (CNV) analysis
- `easyXpress_function.R`
Wrapped function that processes and cleans the raw animal length data with design files using edited _easyXpress_ functions saving figures and files to designated output source 
- `easyXpress_loop.R` 
Unwrapped version of easyXpress_function.R
- `figures.R` 
Generates the main figure for _Genetic differences impact anthelmintic responses across wild Caenorhabditis briggsae strains_ paper
- `orthofinder3.sh` 
Run the _OrthoFinder_ command

