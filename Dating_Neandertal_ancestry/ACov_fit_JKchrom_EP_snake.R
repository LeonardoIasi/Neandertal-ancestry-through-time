
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


grid_fit_function <- function(td_values,tm_values,data.in){
  res_matrix.EP <- matrix(NA,nrow = length(td_values),ncol =  length(tm_values),byrow = T)
  
  try({
    for(i in 1:length(td_values)){
      for(j in 1:length(tm_values)){
        res_fit <- GriD_EP_ACov_amplitude_x_inter_fit_weighted(data.in = data.in,td_fixed = td_values[i],tm_fixed = tm_values[j],A_init = data.in$A_init,C_init = data.in$C_init)
        res_matrix.EP[i,j] <- res_fit$RSS
      }
    }
    colnames(res_matrix.EP) <- tm_values
    rownames(res_matrix.EP) <- td_values
    
    n_DP = length(data.in$y)
    res_matrix_LL.EP <- matrix(NA,nrow = length(td_values),ncol =  length(tm_values),byrow = T)
    
    for(i in 1:length(td_values)){
      for(j in 1:length(tm_values)){
        res_matrix_LL.EP[i,j] <-  - n_DP/2 * log(res_matrix.EP[i,j]) - n_DP/2 * log(2*pi / n_DP) - n_DP/2
        
      }
    }
    
    colnames(res_matrix_LL.EP) <- tm_values
    rownames(res_matrix_LL.EP) <- td_values
    
    return(res_matrix_LL.EP)
  }, silent=F)
  return(res_matrix_LL.EP <- matrix(NA,nrow = 1,ncol =  1,byrow = T))
  
  
}

### input Snakemake

Meta = snakemake@input$meta_data
input_file_folder = snakemake@input$input_folder

map_ = snakemake@wildcards$Map
trunc_ = snakemake@wildcards$truncation
exclude_ = strsplit(snakemake@wildcards$exclude_samples,"-")[[1]]
asc_ = as.character(snakemake@wildcards$Asc)
min_l = as.numeric(snakemake@wildcards$min_l_ACov)
weighted = as.character(snakemake@wildcards$weighted)
gen_time = as.numeric(snakemake@params$gen_time)
step_size_td = as.numeric(snakemake@params$step_size_td)
step_size_tm = as.numeric(snakemake@params$step_size_tm)
max_td = as.numeric(snakemake@params$max_td)
max_tm = as.numeric(snakemake@params$max_tm)
JKchrom = as.numeric(snakemake@wildcards$JKchrom)

output_SP = snakemake@output$output_SP
output_EP = snakemake@output$output_LL_matrix



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

truncation_df = data.frame(sample = c("BK1653","BKBB7240","BKF6620","BKCC7335","OaseNew","UstIshim","ZlatyKunShotgun"),
                           trunc = c(2,5,8,5,25,5,5))

ACov_list <- list()
offset_vector <- c()
sample_names_ACov <- c()
for(i in 1:length(samples)){
  smpl = samples[i]
  if(smpl %in% c("BK1653","BKBB7240","BKF6620","BKCC7335","OaseNew","UstIshim","ZlatyKunShotgun")) {
    print("Sampling longer ACov")
    x <- read.table(paste0(input_file_folder,"longer_distance/outfiles_",smpl,"/",smpl,":",JKchrom))  
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



# SP

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


## EP


print(paste0("calculating JK for EP ACov: ",JKchrom))


# define th grid
td_values_x = unique(c(1,seq(from = step_size_td, to = max_td, by = step_size_td)))
tm_values_x = unique(c(round(max(offset_vector)) + 1, round(max(offset_vector)) + seq(from = step_size_tm, to = max_tm, by = step_size_tm)))



data.in = list(
  dist=lapply(ACov_list, `[[`, 'dist'),
  ACov_res_list=ACov_list,
  y=unlist(lapply(ACov_list, `[[`, 'wcov')),
  weight = unlist(lapply(ACov_list, `[[`, 'weight')),
  offset = offset_vector,
  A_init = rep(0.1,length(offset_vector)),
  C_init = rep(0,length(offset_vector))
)

res_matrix_ACov_fit_LL.EP <- grid_fit_function(td_values_x,tm_values_x,data.in)

res_melt_ACov_fit_LL.EP <- reshape2::melt(res_matrix_ACov_fit_LL.EP)
res_melt_ACov_fit_LL.EP$value0 <- res_melt_ACov_fit_LL.EP$value - max(res_melt_ACov_fit_LL.EP$value)

est_td_EP_ACov_fit_EP=res_melt_ACov_fit_LL.EP$X1[res_melt_ACov_fit_LL.EP$value0 == max(res_melt_ACov_fit_LL.EP$value0)]
est_tm_EP_ACov_fit_EP=res_melt_ACov_fit_LL.EP$X2[res_melt_ACov_fit_LL.EP$value0 == max(res_melt_ACov_fit_LL.EP$value0)]

res_melt_ACov_fit_LL.EP$value0 <- ifelse(res_melt_ACov_fit_LL.EP$value0 < -100,-101,res_melt_ACov_fit_LL.EP$value0)
res_melt_ACov_fit_LL.EP$method <- "ACov_fit"
res_melt_ACov_fit_LL.EP$model <- "ExtendedPulse"
res_melt_ACov_fit_LL.EP$est_tm <- est_tm_EP_ACov_fit_EP
res_melt_ACov_fit_LL.EP$est_td <- est_td_EP_ACov_fit_EP

res_melt_ACov_fit_LL.EP$min_l <- min_l
res_melt_ACov_fit_LL.EP$exclude <- exclude_
res_melt_ACov_fit_LL.EP$map <- map_
res_melt_ACov_fit_LL.EP$weighted <- weighted
res_melt_ACov_fit_LL.EP$asc <- asc_
res_melt_ACov_fit_LL.EP$truncation <- trunc_
res_melt_ACov_fit_LL.EP$gen_time <- gen_time
res_melt_ACov_fit_LL.EP$step_size_tm <- step_size_tm
res_melt_ACov_fit_LL.EP$step_size_td <- step_size_td
res_melt_ACov_fit_LL.EP$max_tm <- max_tm
res_melt_ACov_fit_LL.EP$max_td <- max_td
res_melt_ACov_fit_LL.EP$JKchrom <- JKchrom

write.csv(x = res_melt_ACov_fit_LL.EP,file = as.character(output_EP),quote = F,row.names = F)




