
library("DEoptim")
library(VGAM)
library(tidyverse)
library(DPQ)
library(MASS) 
library(reshape2) 
library(reshape) 
'%notin%' <- Negate('%in%')

# functions

source("~/EMH_Introgression_Project/Introgression_Detection/scripts/Analysis_R_functions.R")


### input Snakemake

Meta = snakemake@input$meta_data

map_ = snakemake@wildcards$Map
trunc_ = snakemake@wildcards$truncation
exclude_ = strsplit(snakemake@wildcards$exclude_samples,"-")[[1]]
asc_ = as.character(snakemake@wildcards$Asc)
min_l = as.numeric(snakemake@wildcards$min_l_ACov)
gen_time = as.numeric(snakemake@params$gen_time)
weighted = as.character(snakemake@wildcards$weighted)


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


input_file_folder = paste0("/mnt/diversity/leonardo_iasi/EMH_Introgression_Project/EMH_dating_analysis/",map_,"_moorjani_et_al/")


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

truncation_df = data.frame(sample = c("BK1653","BKBB7240","BKF6620","BKCC7335","OaseNew","UstIshim","ZlatyKunShotgun"),
                           trunc = c(2,5,8,5,25,5,5))

# SP

ACov_list <- list()
offset_vector <- c()
sample_names_ACov <- c()
for(i in 1:length(samples)){
  smpl = samples[i]
  if(smpl %in% c("BK1653","BKBB7240","BKF6620","BKCC7335","OaseNew","UstIshim","ZlatyKunShotgun")) {
    print("Sampling longer ACov")
    x <- read.table(paste0(input_file_folder,"longer_distance/outfiles_",smpl,"/",smpl))  
  } else{
    x <- read.table(paste0(input_file_folder,"outfiles_",smpl,"/",smpl))
  }
  if(trunc_ == "not_truncated"){
    max_l = 60
  }
  if(trunc_ == "estimated_truncation"){
    max_l = x$V1[which(x$V2 < 0)[1]]
    if(is.na(max_l)){
      max_l = 1
    }
    if( max_l < 1){
      max_l = 1
    }
    max_l = round_any(max_l,1, f = ceiling) 
  }
  if(trunc_ == "fixed_truncation"){
    if(smpl %in% c("BK1653","BKBB7240","BKF6620","BKCC7335","OaseNew","UstIshim","ZlatyKunShotgun")){
      max_l = truncation_df$trunc[truncation_df$sample == smpl]
    } else {
      max_l = 1
    }
  }
  
  x <- x %>% filter(V1 >= min_l, V1 <= max_l) %>% 
    dplyr::rename(dist = V1, wcov = V2) %>% mutate(dist = dist/100)
  
  if(weighted == "yes"){
    x <- x %>% mutate(weight = 1000/length(dist))
  } else {
    x <- x %>% mutate(weight = 1)
  }
  
  time_x <- Joint_Meta$ML_BP_Mean[Joint_Meta$sample_name == smpl] 
  ACov_list[[i]] <- x
  offset_vector <- c(offset_vector,time_x)
  sample_names_ACov <- c(sample_names_ACov,smpl)
}

offset_vector <- offset_vector / gen_time 



data.in = list(
  dist=lapply(ACov_list, `[[`, 'dist'),
  ACov_res_list=ACov_list,
  y=unlist(lapply(ACov_list, `[[`, 'wcov')),
  weight = unlist(lapply(ACov_list, `[[`, 'weight')),
  offset = offset_vector,
  A_init = rep(0.1,length(offset_vector)),
  tm_init = max(offset_vector)+50,
  C_init = rep(0,length(offset_vector))
)


res_SP = Joint_Simple_Pulse_ACov_weighted_fn(data.in = data.in,minFac = 1e-4,maxIter = 5,JointInterceptsFit = F)

tm_CI_SP_ACov <- round(Get_CI_Dating_fn(est = res_SP$tm,n_data = length(data.in$y)))

RSS_SP_ACov <- res_SP$RSS
n_DP = length(data.in$y)
LL_SP_ACov_fit = - n_DP/2 * log(RSS_SP_ACov) - n_DP/2 * log(2*pi / n_DP) - n_DP/2

Sp_res <- data.frame(tm = tm_CI_SP_ACov[1],tm_lower = tm_CI_SP_ACov[2],tm_upper = tm_CI_SP_ACov[3],LL=LL_SP_ACov_fit,
                     method="ACov_fit",model="SimplePulse",min_l = min_l,exclude = exclude_,map=map_,asc=asc_,truncation = trunc_,weighted=weighted,gen_time = gen_time,JKchrom="All")



# SP JKchrom

Sp_res_JK_list <- list()
for(chr in 1:22){
  ACov_list <- list()
  offset_vector <- c()
  sample_names_ACov <- c()
  for(i in 1:length(samples)){
    smpl = samples[i]
    if(smpl %in% c("BK1653","BKBB7240","BKF6620","BKCC7335","OaseNew","UstIshim","ZlatyKunShotgun")) {
      print("Sampling longer ACov")
      x <- read.table(paste0(input_file_folder,"longer_distance/outfiles_",smpl,"/",smpl,":",chr))  
    } else{
      x <- read.table(paste0(input_file_folder,"outfiles_",smpl,"/",smpl))
    }
    if(trunc_ == "not_truncated"){
      max_l = 60
    }
    if(trunc_ == "estimated_truncation"){
      max_l = x$V1[which(x$V2 < 0)[1]]
      if(is.na(max_l)){
        max_l = 1
      }
      if( max_l < 1){
        max_l = 1
      }
      max_l = round_any(max_l,1, f = ceiling) 
    }
    if(trunc_ == "fixed_truncation"){
      if(smpl %in% c("BK1653","BKBB7240","BKF6620","BKCC7335","OaseNew","UstIshim","ZlatyKunShotgun")){
        max_l = truncation_df$trunc[truncation_df$sample == smpl]
      } else {
        max_l = 1
      }
    }
    
    x <- x %>% filter(V1 >= min_l, V1 <= max_l) %>% 
      dplyr::rename(dist = V1, wcov = V2) %>% mutate(dist = dist/100)
    
    if(weighted == "yes"){
      x <- x %>% mutate(weight = 1000/length(dist))
    } else {
      x <- x %>% mutate(weight = 1)
    }
    
    time_x <- Joint_Meta$ML_BP_Mean[Joint_Meta$sample_name == smpl] 
    ACov_list[[i]] <- x
    offset_vector <- c(offset_vector,time_x)
    sample_names_ACov <- c(sample_names_ACov,smpl)
  }
  
  offset_vector <- offset_vector / gen_time 
  
  
  
  data.in = list(
    dist=lapply(ACov_list, `[[`, 'dist'),
    ACov_res_list=ACov_list,
    y=unlist(lapply(ACov_list, `[[`, 'wcov')),
    weight = unlist(lapply(ACov_list, `[[`, 'weight')),
    offset = offset_vector,
    A_init = rep(0.1,length(offset_vector)),
    tm_init = max(offset_vector)+50,
    C_init = rep(0,length(offset_vector))
  )
  
  
  res_SP = Joint_Simple_Pulse_ACov_weighted_fn(data.in = data.in,minFac = 1e-4,maxIter = 5,JointInterceptsFit = F)
  
  tm_CI_SP_ACov <- round(Get_CI_Dating_fn(est = res_SP$tm,n_data = length(data.in$y)))
  
  RSS_SP_ACov <- res_SP$RSS
  n_DP = length(data.in$y)
  LL_SP_ACov_fit = - n_DP/2 * log(RSS_SP_ACov) - n_DP/2 * log(2*pi / n_DP) - n_DP/2
  
  Sp_res_JK_list[[chr]] <- data.frame(tm = tm_CI_SP_ACov[1],tm_lower = tm_CI_SP_ACov[2],tm_upper = tm_CI_SP_ACov[3],LL=LL_SP_ACov_fit,
                       method="ACov_fit",model="SimplePulse",min_l = min_l,exclude = exclude_,map=map_,asc=asc_,truncation = trunc_,weighted=weighted,gen_time = gen_time,JKchrom=chr)
  
  
  
}

Sp_res_JK <- do.call("rbind",Sp_res_JK_list)

Sp_res_all <- rbind(Sp_res,Sp_res_JK)

write.csv(x = Sp_res_all,file = as.character(output_SP),quote = F,row.names = F)

