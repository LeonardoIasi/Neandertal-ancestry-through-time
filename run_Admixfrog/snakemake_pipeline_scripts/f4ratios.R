library(tidyverse)
library(plyr)
library(yaml)
setwd("~/EMH_Introgression_Project/Introgression_Detection/scripts/plotting/")

panel = snakemake@wildcards$panel
merged_alleles = snakemake@input$merged_all_alleles
deam = as.logical(snakemake@params$deam)
calc_DEN = as.logical(snakemake@params$calc_DEN)
output = snakemake@output$f4r_all_alleles_wide


f4r <- function(df, X, A, B, C, D){
  X <- enquo(X)
  A <- enquo(A)
  B <- enquo(B)
  C <- enquo(C)
  D <- enquo(D)


  df %>% dplyr::select(chrom, pos, (!!X), (!!A), (!!B), (!!C), (!!D) ) %>%
    dplyr::rename(
      X = (!!X),
      A = (!!A),
      B = (!!B),
      C = (!!C),
      D = (!!D)
    ) %>%
    drop_na() %>%
    summarize(f4r = sum((A-D)*(X-C)) / sum((A-D) * (B-C))) %>% return()
}

NEA_f4r <- function(sample_name,data,filter_distance=NULL,DEN_filter=F){
  print(paste("Processing: ",sample_name,sep = ""))
  ref_sample = data  %>%   dplyr::rename(X = sample_name) %>%
    dplyr::select(chrom, pos, ALT, VIN, DEN, AFR, PAN, X  ) %>% drop_na()

  if(!is.null(filter_distance)){
    ref_sample %>% group_by(chrom) %>% dplyr::summarise(distance = abs(pos - c(pos[2:length(pos)],1000))) %>% drop_na() -> distances_between_snps
    ref_sample$distances_between_snps <-  distances_between_snps$distance
    ref_sample = ref_sample %>% filter(distances_between_snps >= filter_distance)
  }

  if(DEN_filter ==T){
    ref_sample = ref_sample %>% filter(DEN == 0)
  }

  x_f4 = f4r(df=ref_sample,X = "X",A = "ALT",B = "VIN",C = "AFR",D = "PAN")
  x_f4 = cbind(sample_name,x_f4)
  return(x_f4)
}

DEN_f4r <- function(sample_name,data,filter_distance=NULL){
  print(paste("Processing: ",sample_name,sep = ""))
  ref_sample = data  %>%   dplyr::rename(X = sample_name) %>%
    dplyr::select(chrom, pos, EAS, VIN, DEN, AFR, X  ) %>% drop_na()

  if(!is.null(filter_distance)){
    ref_sample %>% group_by(chrom) %>% dplyr::summarise(distance = abs(pos - c(pos[2:length(pos)],1000))) %>% drop_na() -> distances_between_snps
    ref_sample$distances_between_snps <-  distances_between_snps$distance
    ref_sample = ref_sample %>% filter(distances_between_snps >= filter_distance)
  }

  x_f4 = f4r(df=ref_sample,X = "X",A = "VIN",B = "DEN",C = "EAS",D = "AFR")
  x_f4 = cbind(sample_name,x_f4)
  return(x_f4)
}

All_AF_input_all <- read_csv(merged_alleles)

sample_names <-   yaml::read_yaml("~/EMH_Introgression_Project/Introgression_Detection/config/panels_EMH_Introgression_Project.yaml")

sample_names = sample_names[["panels"]][[panel]]
#sample_names = gsub("-",".",sample_names)

All_AF_input_NEA = as.data.frame(matrix(unlist(lapply(sample_names,function(x) NEA_f4r(x,All_AF_input_all,filter_distance=NULL,DEN_filter=F))),ncol = 2,byrow = T))
colnames(All_AF_input_NEA) <- c("sample","f4r")
All_AF_input_NEA$data <- "all"
All_AF_input_NEA$Ancestry <- "NEA"

if(calc_DEN == T){
  All_AF_input_DEN = as.data.frame(matrix(unlist(lapply(sample_names,function(x) DEN_f4r(x,All_AF_input_all,filter_distance=NULL))),ncol = 2,byrow = T))
  colnames(All_AF_input_DEN) <- c("sample","f4r")
  All_AF_input_DEN$data <- "all"
  All_AF_input_DEN$Ancestry <- "DEN"
  if(deam == T){
    merged_alleles = snakemake@input$merged_deam_only
    All_AF_input_deam <- read.csv(merged_alleles)

    All_AF_input_NEA_deam = as.data.frame(matrix(unlist(lapply(sample_names,function(x) NEA_f4r(x,All_AF_input_deam,filter_distance=NULL,DEN_filter=F))),ncol = 2,byrow = T))
    colnames(All_AF_input_NEA_deam) <- c("sample","f4r")
    All_AF_input_NEA_deam$data <- "deam"
    All_AF_input_NEA_deam$Ancestry <- "NEA"

    All_AF_input_DEN_deam = as.data.frame(matrix(unlist(lapply(sample_names,function(x) DEN_f4r(x,All_AF_input_deam,filter_distance=NULL))),ncol = 2,byrow = T))
    colnames(All_AF_input_DEN_deam) <- c("sample","f4r")
    All_AF_input_DEN_deam$data <- "deam"
    All_AF_input_DEN_deam$Ancestry <- "DEN"

    Ancestry_estimates = plyr::join_all(list(All_AF_input_NEA,All_AF_input_DEN,All_AF_input_NEA_deam,All_AF_input_DEN_deam), by='sample', type='inner')
    Ancestry_estimates_long = rbind(All_AF_input_NEA,All_AF_input_DEN,All_AF_input_NEA_deam,All_AF_input_DEN_deam)
    write.csv(x = Ancestry_estimates,file = output,quote = F,row.names = F)
    write.csv(x = Ancestry_estimates_long,file = paste(substr(output,1,nchar(output)-4),"_long.csv",sep=""),quote = F,row.names = F)

  }else{
    Ancestry_estimates = plyr::join_all(list(All_AF_input_NEA,All_AF_input_DEN), by='sample', type='inner')
    Ancestry_estimates_long = rbind(All_AF_input_NEA,All_AF_input_DEN)
    write.csv(x = Ancestry_estimates,file = output,quote = F,row.names = F)
    write.csv(x = Ancestry_estimates_long,file = paste(substr(output,1,nchar(output)-4),"_long.csv",sep=""),quote = F,row.names = F)
  }
}else{
  if(deam == T){
    merged_alleles = snakemake@input$merged_deam_only
    All_AF_input_deam <- read.csv(merged_alleles)

    All_AF_input_NEA_deam = as.data.frame(matrix(unlist(lapply(sample_names,function(x) NEA_f4r(x,All_AF_input_deam,filter_distance=NULL,DEN_filter=F))),ncol = 2,byrow = T))
    colnames(All_AF_input_NEA_deam) <- c("sample","f4r")
    All_AF_input_NEA_deam$data <- "deam"
    All_AF_input_NEA_deam$Ancestry <- "NEA"


    Ancestry_estimates = plyr::join_all(list(All_AF_input_NEA,All_AF_input_NEA_deam), by='sample', type='inner')
    Ancestry_estimates_long = rbind(All_AF_input_NEA,All_AF_input_NEA_deam)
    write.csv(x = Ancestry_estimates,file = output,quote = F,row.names = F)
    write.csv(x = Ancestry_estimates_long,file = paste(substr(output,1,nchar(output)-4),"_long.csv",sep=""),quote = F,row.names = F)

  }else{
    Ancestry_estimates = plyr::join_all(list(All_AF_input_NEA), by='sample', type='inner')
    Ancestry_estimates_long = rbind(All_AF_input_NEA)
    write.csv(x = Ancestry_estimates,file = output,quote = F,row.names = F)
    write.csv(x = Ancestry_estimates_long,file = paste(substr(output,1,nchar(output)-4),"_long.csv",sep=""),quote = F,row.names = F)
  }

}
