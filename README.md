# Genomic-Dominance-Relationship-Matrix
R scripts for calculating the centered and normalised dominance relationship matrices. The scripts are expecting the data to be in {-1, 0, 1} notation and will convert it to {0, 1, 2} notation for the dominance relationship matrices and {0,0.5,1} notation for the Yang method. The program will work if data in the wrong notation are supplied, but it will do so without a warning, and the results are not standard methods in quantitative genetic analysis.

Note that missing data must be imputed. Only the Van Raden method will impute data based on the marker mode. 

Details on methodologies can be found in: 

Alireza Nazarian and Salvador Alejandro Gezan. 2016. GenoMatrix: a Software Package for Pedigree-based and Genomic Prediction Analyses on Complex Traits. J Hered. doi:10.1093/jhered/esw020
