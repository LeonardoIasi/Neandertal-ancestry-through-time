library(tidyverse)
library(plyr)
setwd("~/EMH_Introgression_Project/Introgression_Detection/scripts/")

ref_file = snakemake@input$ref

input_all = snakemake@input$all_alleles
deam = as.logical(snakemake@params$deam)
calc_DEN = as.logical(snakemake@params$calc_DEN)
output = snakemake@output$merged_all_alleles

if(calc_DEN == T){
  ref = read_csv(ref_file,col_types = list(chrom=col_character()))
  ref = ref %>% 
    mutate(ALT = ALT_alt/(ALT_alt+ALT_ref)) %>%
    mutate(VIN = VIN_alt/(VIN_alt+VIN_ref)) %>%
    mutate(DEN = DEN_alt/(DEN_alt+DEN_ref)) %>%
    mutate(AFR = AFR_alt/(AFR_alt+AFR_ref)) %>%
    mutate(EAS = EAS_alt/(EAS_alt+EAS_ref)) %>%
    mutate(EUR = EUR_alt/(EUR_alt+EUR_ref)) %>%
    mutate(PAN = PAN_alt/(PAN_alt+PAN_ref)) %>%
    dplyr::select(chrom, pos, ALT, VIN, DEN, AFR, EAS, EUR, PAN) %>% mutate_all(function(x) ifelse(is.nan(x), NA, x))
  
  
  All_AF_input_all =  lapply(input_all, function(x) read_csv(x,col_types = list(chrom=col_character()))) %>% purrr::reduce(full_join, by = c("chrom","pos")) %>% full_join(ref,.,by = c("chrom","pos"))
  write.csv(x=All_AF_input_all,file = output,quote = F,row.names = F)
  if(deam ==T){
    input_all_deam = snakemake@input$deam_only
    All_AF_input_deam = lapply(input_all_deam, function(x) read_csv(x,col_types = list(chrom=col_character()))) %>% purrr::reduce(full_join, by = c("chrom","pos")) %>% full_join(ref,.,by = c("chrom","pos"))
    write.csv(x=All_AF_input_deam,file = paste(substr(output,1,nchar(output)-4),"_deam.csv",sep=""),quote = F,row.names = F)
  }
}else{
  ref = read_csv(ref_file,col_types = list(chrom=col_character())) 
  ref = ref %>% 
    mutate(ALT = ALT_alt/(ALT_alt+ALT_ref)) %>%
    mutate(VIN = VIN_alt/(VIN_alt+VIN_ref)) %>%
    mutate(DEN = DEN_alt/(DEN_alt+DEN_ref)) %>%
    mutate(AFR = AFR_alt/(AFR_alt+AFR_ref)) %>%
    mutate(PAN = PAN_alt/(PAN_alt+PAN_ref)) %>%
    dplyr::select(chrom, pos, ALT, VIN, DEN, AFR, PAN) %>% mutate_all(function(x) ifelse(is.nan(x), NA, x))

  
  All_AF_input_all =  lapply(input_all, function(x) read_csv(x,col_types = list(chrom=col_character()))) %>% purrr::reduce(full_join, by = c("chrom","pos")) %>% full_join(ref,.,by = c("chrom","pos"))
  write.csv(x=All_AF_input_all,file = output,quote = F,row.names = F)
  if(deam ==T){
    input_all_deam = snakemake@input$deam_only
    All_AF_input_deam = lapply(input_all_deam, function(x) read_csv(x,col_types = list(chrom=col_character()))) %>% purrr::reduce(full_join, by = c("chrom","pos")) %>% full_join(ref,.,by = c("chrom","pos"))
    write.csv(x=All_AF_input_deam,file = paste(substr(output,1,nchar(output)-4),"_deam.csv",sep=""),quote = F,row.names = F)
  }
  
}


