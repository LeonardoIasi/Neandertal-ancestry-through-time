# Neandertal-ancestry-through-time
This repository contains the script used in the Neandertal ancestry through time project. All code that was used in the project is uploaded as is, if not stated otherwise.

## Pipline for calling Neandertal ancestry
In the [run_Admixfrog](https://github.com/LeonardoIasi/Neandertal-ancestry-through-time/edit/main/run_Admixfrog) folder you can find the snakemake pipline that were used to run [admixfrog](https://github.com/BenjaminPeter/admixfrog) (Peter 2021) on the genomes to call Neandertal ancestry, generat summary statistics and match the segments to the reference seuqunces. 

## Dating Neandertal ancestry 
In the [Dating_Neandertal_ancestry](https://github.com/LeonardoIasi/Neandertal-ancestry-through-time/tree/main/Dating_Neandertal_ancestry) folder you can find the scripts that were used for the dating of the Neandertal gene flow into modern humans using individuals older than 20 ky. Most scripts are based on the ancestry covariance curves from the [Probabilistic-Neanderthal-admixture-dating](https://github.com/MoorjaniLab/Probabilistic-Neanderthal-admixture-dating/tree/master), a method by Moorjani et al. 2016. 

## Analysis
In the [Analysis](https://github.com/LeonardoIasi/Neandertal-ancestry-through-time/tree/main/Analysis) folder you can find the R markdown scripts that were used for the analysis of the called segments from admixfrog. These files are uploaded with the original code however, they are made to take in the files one can download from the associated Dryad repo to make the rerun of the analysis easier. The QC files are uploaded as is.
