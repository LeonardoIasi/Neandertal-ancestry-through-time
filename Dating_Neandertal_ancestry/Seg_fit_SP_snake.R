
library("DEoptim")
library(VGAM)
library(tidyverse)
library(DPQ)
library(MASS) 
library(reshape2) 
library(reshape) 
'%notin%' <- Negate('%in%')

# functions

source("Neandertal-ancestry-through-time/Analysis/Analysis_R_functions.R")


### input Snakemake

Meta = snakemake@input$meta_data
input_file_folder = snakemake@input$input_folder

map_ = snakemake@wildcards$Map
exclude_ = strsplit(snakemake@wildcards$exclude_samples,"-")[[1]]
asc_ = as.character(snakemake@wildcards$Asc)
min_l = as.numeric(snakemake@wildcards$min_l_Seg)
gen_time = as.numeric(snakemake@params$gen_time)

output_SP = snakemake@output$output_SP


# Meta
Joint_Meta <- read.csv(Meta) %>% mutate(PCA_pops = ifelse(ML_BP_Mean == 0 ,pop,sample_name),Cov = ifelse(ML_BP_Mean == 0 ,ave_cov_reported,ave_cov_array),DataType = ifelse(Data_source == "Shotgun" | Data_source == "Genotypes","Shotgun","Captured" )) %>% filter(SGDP_Superpop != "Africa" )

min_age_BP = 20000

# dates from Djakovic et al 2022 for Spain and France
Protoaurignacian_start = 42653
Chatelperronian_extinction = 39798
Neandertal_extinction_Djakovic = 40457

# Higham et al 2014
TransitionalIndustry_start = 43000
Nenadertal_extinction_Higham = 39220

Archeaology_first_humans = max(c(Protoaurignacian_start,TransitionalIndustry_start)) / gen_time
Archeaology_last_neandertals = min(c(Chatelperronian_extinction,Neandertal_extinction_Djakovic,Nenadertal_extinction_Higham)) / gen_time

samples = Joint_Meta$sample_name[Joint_Meta$ML_BP_Mean >= min_age_BP]

if(exclude_ == "None"){
  samples_data = Joint_Meta %>% filter(ML_BP_Mean >= min_age_BP) %>% 
    dplyr::select(sample_name,ML_BP_Mean,ML_BP_Lower,ML_BP_Higher,DataType) 
  exclude_ = ""
  print("No sample excluded")
} else{
  samples_data = Joint_Meta %>% filter(ML_BP_Mean >= min_age_BP) %>% 
    dplyr::select(sample_name,ML_BP_Mean,ML_BP_Lower,ML_BP_Higher,DataType) %>% filter(sample_name %notin% exclude_)
  print(paste("Samples :",exclude_,"excluded"))
  exclude_ = paste(exclude_,collapse  = "_")
}

if(asc_ == "All"){
  samples_data = samples_data
}
if(asc_ == "Shotgun"){
  samples_data = samples_data %>% filter(DataType == "Shotgun")
}
if(asc_ == "Captured"){
  samples_data = samples_data %>% filter(DataType == "Captured")
}

samples <- samples_data$sample_name

dir_path_EMH <- input_file_folder

chrom_ <- c("1" ,"2" , "3" , "4",  "5" , "6"  ,"7" , "8" , "9" , "10" ,"11", "12" ,"13","14" ,"15", "16" ,"17" ,"18", "19" , "20" ,"21", "22"  )
states = c("NEA")

EMH_rle_all <- read.csv(input_file_folder) %>% filter(genetic_map == !!map, target == !!states, type == "state", chrom %in% !!chrom_) %>%
  drop_na() %>% filter(map_len >= 0.05) %>% inner_join(.,Joint_Meta[,c("sample_name","ML_BP_Mean","ML_BP_Higher","ML_BP_Lower")],by=c("sample"="sample_name"))


# SP all

data.in_seg = EMH_rle_all %>% filter(ML_BP_Mean >= min_age_BP) %>% 
  filter(map_len >= min_l) %>% 
  filter(sample %in% samples)

Segment_fit_Data <- list(
  seg_len = as.numeric(data.in_seg$map_len),
  ind = as.numeric(as.factor(data.in_seg$sample)),
  sample_time = as.numeric(data.in_seg$ML_BP_Mean)/gen_time
)

l = Segment_fit_Data$seg_len
offset_vector_segment = Segment_fit_Data$sample_time
lower_trunc = min_l
upper_trunc = max(l)
sample_names_Segments <- levels(as.factor(data.in_seg$sample))


tm_lower = max(Segment_fit_Data$sample_time) + 1
tm_upper = max(Segment_fit_Data$sample_time) + 1000


res = optim(c(tm_lower + 10), fnx_SP_seg, method="L-BFGS-B",lower = c(tm_lower),upper = c(tm_upper))
par =res$par
tm_SP_seg_fit <- round(Get_CI_Dating_fn(est = par,n_data = length(l)))
LL_SP_seg_fit = res$value

Sp_res <- data.frame(tm = tm_SP_seg_fit[1],tm_lower = tm_SP_seg_fit[2],tm_upper = tm_SP_seg_fit[3],LL=LL_SP_seg_fit,
                     method="Segments",model="SimplePulse",min_l = min_l,exclude = exclude_,map=map_,asc=asc_,truncation = NA,weighted=NA,gen_time = gen_time,JKchrom="All")

# SP JKchrom


Sp_res_JK_list <- list()
for(chr in 1:22){
  data.in_seg = EMH_rle_all %>% filter(ML_BP_Mean >= min_age_BP) %>% 
    filter(map_len >= min_l) %>% 
    filter(sample %in% samples) %>%
    filter(chrom != chr)
  
  Segment_fit_Data <- list(
    seg_len = as.numeric(data.in_seg$map_len),
    ind = as.numeric(as.factor(data.in_seg$sample)),
    sample_time = as.numeric(data.in_seg$ML_BP_Mean)/gen_time
  )
  
  l = Segment_fit_Data$seg_len
  offset_vector_segment = Segment_fit_Data$sample_time
  lower_trunc = min_l
  upper_trunc = max(l)
  sample_names_Segments <- levels(as.factor(data.in_seg$sample))
  
  
  tm_lower = max(Segment_fit_Data$sample_time) + 1
  tm_upper = max(Segment_fit_Data$sample_time) + 1000
  
  
  res = optim(c(tm_lower + 10), fnx_SP_seg, method="L-BFGS-B",lower = c(tm_lower),upper = c(tm_upper))
  par =res$par
  tm_SP_seg_fit <- round(Get_CI_Dating_fn(est = par,n_data = length(l)))
  LL_SP_seg_fit = res$value
  
  Sp_res_JK_list[[chr]] <- data.frame(tm = tm_SP_seg_fit[1],tm_lower = tm_SP_seg_fit[2],tm_upper = tm_SP_seg_fit[3],LL=LL_SP_seg_fit,
                       method="Segments",model="SimplePulse",min_l = min_l,exclude = exclude_,map=map_,asc=asc_,truncation = NA,weighted=NA,gen_time = gen_time,JKchrom=chr)
  

  
  
}

Sp_res_JK <- do.call("rbind",Sp_res_JK_list)

Sp_res_all <- rbind(Sp_res,Sp_res_JK)

write.csv(x = Sp_res_all,file = as.character(output_SP),quote = F,row.names = F)
