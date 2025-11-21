# Neandertal-ancestry-through-time
This repository contains the script used in the Neandertal ancestry through time project. All code that was used in the project is uploaded as is, if not stated otherwise.
Most input data for the code can be found on [Dryad](https://datadryad.org/dataset/doi:10.5061/dryad.zw3r228gg).

## Pipline for calling Neandertal ancestry
In the [run_Admixfrog](https://github.com/LeonardoIasi/Neandertal-ancestry-through-time/edit/main/run_Admixfrog) folder you can find the snakemake pipline that were used to run [admixfrog](https://github.com/BenjaminPeter/admixfrog) (Peter 2021) on the genomes to call Neandertal ancestry, generat summary statistics and match the segments to the reference seuqunces. 

## Dating Neandertal ancestry 
To get the Ancestry Covariance (method based on Moorjani et al. 2016) and singe sample admixture time follow the description in the [Dating-and_Outlier-Scan](https://github.com/LeonardoIasi/Dating-and-Outlier-Scan/tree/d7c23214e2e9e06a16a6e3f94cdb3416bc76f82f), subfolder Genetic_Dating. The input requires ancestry informative (ascertained) SNP information in EIGEN format and a genetic map. You can find the Covariance date for this study on the dryad repositiry [Dryad](https://datadryad.org/dataset/doi:10.5061/dryad.zw3r228gg) in the Ancestry_Covariance_Data_Shared_map.zip folder.
In the [Dating_Neandertal_ancestry](https://github.com/LeonardoIasi/Neandertal-ancestry-through-time/tree/main/Dating_Neandertal_ancestry) folder you can find the scripts that use the Ancestry Covariance for jointly dating the Neandertal gene flow into modern humans using individuals older than 20 ky. 

## Analysis
In the [Analysis](https://github.com/LeonardoIasi/Neandertal-ancestry-through-time/tree/main/Analysis) folder you can find the R markdown scripts that were used for the analysis of the called segments from admixfrog. These files are uploaded with the original code however, they are made to take in the files one can download from the associated Dryad repo to make the rerun of the analysis easier. The QC files are uploaded as is.

For any questions plaese send a e-mail to: leonardo_iasi@eva.mpg.de
