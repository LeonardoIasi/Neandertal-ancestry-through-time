
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

sample_set = strsplit(snakemake@wildcards$sample_set,"-")[[1]]
map_ = snakemake@wildcards$Map
trunc_ = snakemake@wildcards$truncation
min_l = as.numeric(snakemake@wildcards$min_l_ACov)
gen_time = as.numeric(snakemake@params$gen_time)
step_size_td = as.numeric(snakemake@params$step_size_td)
step_size_tm = as.numeric(snakemake@params$step_size_tm)
max_td = as.numeric(snakemake@params$max_td)
max_tm = as.numeric(snakemake@params$max_tm)
asc_ = "All"
weighted = "no"
JKchrom = as.numeric(snakemake@wildcards$JKchrom)

output_SP = snakemake@output$output_SP
output_EP = snakemake@output$output_EP
output_EP_grid = snakemake@output$output_EP_grid

output_summary <- snakemake@output$output_summary

print(trunc_)

# Meta
Joint_Meta <- read.csv(Meta) %>% mutate(PCA_pops = ifelse(ML_BP_Mean == 0 ,pop,sample_name),Cov = ifelse(ML_BP_Mean == 0 ,ave_cov_reported,ave_cov_array),DataType = ifelse(Data_source == "Shotgun" | Data_source == "Genotypes","Shotgun","Captured" )) %>% filter(SGDP_Superpop != "Africa" )

min_age_BP = 20000


samples = Joint_Meta$sample_name[Joint_Meta$ML_BP_Mean >= min_age_BP]


input_file_folder = paste0("/mnt/diversity/leonardo_iasi/EMH_Introgression_Project/EMH_dating_analysis/",map_,"_moorjani_et_al/")



samples_data = Joint_Meta %>% filter(ML_BP_Mean >= min_age_BP) %>% 
  dplyr::select(sample_name,ML_BP_Mean,ML_BP_Lower,ML_BP_Higher,DataType) 



samples <- samples_data$sample_name[samples_data$sample_name %in% sample_set]

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
  
  
  time_x <- Joint_Meta$ML_BP_Mean[Joint_Meta$sample_name == smpl] 
  ACov_list[[i]] <- x
  offset_vector <- c(offset_vector,time_x)
  sample_names_ACov <- c(sample_names_ACov,smpl)
}

offset_vector <- offset_vector / gen_time 

# SP

Single_Sample_SP_ACov <- data.frame(sample = rep(NA,length(unique(sample_names_ACov))),
                                    tm_SP = rep(NA,length(unique(sample_names_ACov))),
                                    A_SP = rep(NA,length(unique(sample_names_ACov))),
                                    c_SP = rep(NA,length(unique(sample_names_ACov))),
                                    n_SP = rep(NA,length(unique(sample_names_ACov))),
                                    LL_SP = rep(NA,length(unique(sample_names_ACov))))


for(i in 1:length(unique(sample_names_ACov))){
  ind = unique(sample_names_ACov)[i]
  
  dist=ACov_list[[i]]$dist
  y=ACov_list[[i]]$wcov
  offset = 0
  A_init = 0.1
  C_init = 0.1
  tm_init = 10
  
  res_fit <- nls(y ~ fnx_SP_singlefit(A,tm,C,0,dist), start=list(A=A_init,tm=tm_init,C=C_init), algorithm="port",
                 lower=c(0,1,-1), upper=c(1,5000,1),
                 control=list(maxiter=50,
                              warnOnly=TRUE,minFactor=0.0004))
  
  RSS = sum((y - fnx_SP_singlefit(coef(res_fit)[1],coef(res_fit)[2],coef(res_fit)[3],0,dist))^2)
  n_DP = length(dist)
  LL = - n_DP/2 * log(RSS) - n_DP/2 * log(2*pi / n_DP) - n_DP/2
  Single_Sample_SP_ACov[i,1] = ind
  Single_Sample_SP_ACov[i,2] = coef(res_fit)[2]
  Single_Sample_SP_ACov[i,3] = coef(res_fit)[1]
  Single_Sample_SP_ACov[i,4] = coef(res_fit)[3]
  Single_Sample_SP_ACov[i,5] = length(dist)
  Single_Sample_SP_ACov[i,6] = LL
  
  
}

Single_Sample_SP_ACov$JKchrom <- JKchrom

# EP fit with nls

Single_Sample_EP_ACov <- data.frame(sample = rep(NA,length(unique(sample_names_ACov))),
                                    tm_EP = rep(NA,length(unique(sample_names_ACov))),
                                    td_EP = rep(NA,length(unique(sample_names_ACov))),
                                    A_EP = rep(NA,length(unique(sample_names_ACov))),
                                    c_EP = rep(NA,length(unique(sample_names_ACov))),
                                    n_EP = rep(NA,length(unique(sample_names_ACov))),
                                    LL_EP = rep(NA,length(unique(sample_names_ACov))),
                                    converged = rep(NA,length(unique(sample_names_ACov))))


for(i in 1:length(unique(sample_names_ACov))){
  ind = unique(sample_names_ACov)[i]
  
  dist=ACov_list[[i]]$dist
  y=ACov_list[[i]]$wcov
  offset = 0
  A_init = 0.1
  C_init = 0.1
  tm_init = 50
  td_init = 5
  
  
  res_fit <- nls(y ~ fnx_EP_singlefit(A,td,tm,C,0,dist), start=list(A=A_init,td=td_init,tm=tm_init,C=C_init), algorithm="port",
                 lower=c(0,1,1,-1), upper=c(1,500,1000,1),
                 control=list(maxiter=100,
                              warnOnly=T,minFactor=0.0004))
  
  RSS = sum((y - fnx_EP_singlefit(coef(res_fit)[1],coef(res_fit)[2],coef(res_fit)[3],coef(res_fit)[4],0,dist))^2)
  n_DP = length(dist)
  LL = - n_DP/2 * log(RSS) - n_DP/2 * log(2*pi / n_DP) - n_DP/2
  Single_Sample_EP_ACov[i,1] = ind
  Single_Sample_EP_ACov[i,2] = coef(res_fit)[3]
  Single_Sample_EP_ACov[i,3] = coef(res_fit)[2]
  Single_Sample_EP_ACov[i,4] = coef(res_fit)[1]
  Single_Sample_EP_ACov[i,5] = coef(res_fit)[4]
  Single_Sample_EP_ACov[i,6] = length(dist)
  Single_Sample_EP_ACov[i,7] = LL
  Single_Sample_EP_ACov[i,8] = res_fit$convInfo$isConv
  
  
}

Single_Sample_EP_ACov$JKchrom <- JKchrom

# EP fit with grid

EP_fit_list <- list()
for(x in 1:length(unique(sample_names_ACov))){
  ind = unique(sample_names_ACov)[x]
  print(ind)
  # define th grid
  td_values = unique(c(1,seq(from = step_size_td, to = max_td, by = step_size_td)))
  tm_values = unique(c(round(max(0)) + 1, round(max(0)) + seq(from = step_size_td, to = ,max_tm, by = step_size_td)))
  
  
  
  data.in = list(
    dist=ACov_list[[x]]$dist,
    ACov_res_list=list(ACov_list[[x]]),
    y=ACov_list[[x]]$wcov,
    weight = 1,
    offset = 0,
    A_init = 0.1,
    C_init = 0.1
  )
  
  res_matrix_ACov_fit.EP <- matrix(NA,nrow = length(td_values),ncol =  length(tm_values),byrow = T)
  res_amplitude.EP <- matrix(NA,nrow = length(td_values) * length(tm_values),ncol =  length(offset_vector),byrow = T)
  res_C.EP <- matrix(NA,nrow = length(td_values) * length(tm_values),ncol =  length(offset_vector),byrow = T)
  
  
  
  
  for(i in 1:length(td_values)){
    for(j in 1:length(tm_values)){
      y = data.in$y
      fitx <- nls(y ~ fnx_EP_singlefit(A,td=td_values[i],tm=tm_values[j],C,offset=0,dist=data.in$dist,td_invers = F),start=list(A=data.in$A_init,C=data.in$C_init), algorithm="port",
                  lower=c(1e-6,0), upper=c(1,1),
                  control=list(warnOnly=TRUE))
      RSS <- sum(((data.in$y - fnx_EP_singlefit(A=coef(fitx)[[1]],td=td_values[i],tm=tm_values[j],C=coef(fitx)[[2]],offset=0,dist=data.in$dist,td_invers = F))^2))
      res_matrix_ACov_fit.EP[i,j] <- RSS
    }
  }
  
  colnames(res_matrix_ACov_fit.EP) <- tm_values
  rownames(res_matrix_ACov_fit.EP) <- td_values
  
  n_DP = length(data.in$y)
  res_matrix_ACov_fit_LL.EP <- matrix(NA,nrow = length(td_values),ncol =  length(tm_values),byrow = T)
  
  for(i in 1:length(td_values)){
    for(j in 1:length(tm_values)){
      res_matrix_ACov_fit_LL.EP[i,j] <-  - n_DP/2 * log(res_matrix_ACov_fit.EP[i,j]) - n_DP/2 * log(2*pi / n_DP) - n_DP/2
      
    }
  }
  
  colnames(res_matrix_ACov_fit_LL.EP) <- tm_values
  rownames(res_matrix_ACov_fit_LL.EP) <- td_values
  
  
  res_all_normal_fits <- list(res_matrix_ACov_fit.EP,res_matrix_ACov_fit_LL.EP)
  
  
  
  # ACov fit
  res_matrix_ACov_fit_LL.EP = res_all_normal_fits[[2]]
  res_melt_ACov_fit_LL.EP <- reshape2::melt(res_matrix_ACov_fit_LL.EP)
  res_melt_ACov_fit_LL.EP$value0 <- res_melt_ACov_fit_LL.EP$value - max(res_melt_ACov_fit_LL.EP$value)
  
  est_td_EP_ACov_fit_EP=res_melt_ACov_fit_LL.EP$Var1[res_melt_ACov_fit_LL.EP$value0 == max(res_melt_ACov_fit_LL.EP$value0)]
  est_tm_EP_ACov_fit_EP=res_melt_ACov_fit_LL.EP$Var2[res_melt_ACov_fit_LL.EP$value0 == max(res_melt_ACov_fit_LL.EP$value0)]
  
  res_melt_ACov_fit_LL.EP$value0 <- ifelse(res_melt_ACov_fit_LL.EP$value0 < -100,-101,res_melt_ACov_fit_LL.EP$value0)
  res_melt_ACov_fit_LL.EP$method <- "ACov_fit"
  res_melt_ACov_fit_LL.EP$model <- "ExtendedPulse"
  res_melt_ACov_fit_LL.EP$est_tm <- est_tm_EP_ACov_fit_EP
  res_melt_ACov_fit_LL.EP$est_td <- est_td_EP_ACov_fit_EP
  
  res_melt_ACov_fit_LL.EP$min_l <- min_l
  res_melt_ACov_fit_LL.EP$map <- map_
  res_melt_ACov_fit_LL.EP$weighted <- weighted
  res_melt_ACov_fit_LL.EP$asc <- asc_
  res_melt_ACov_fit_LL.EP$truncation <- trunc_
  res_melt_ACov_fit_LL.EP$gen_time <- gen_time
  res_melt_ACov_fit_LL.EP$ind <- ind
  res_melt_ACov_fit_LL.EP$max_td <- max_td
  res_melt_ACov_fit_LL.EP$max_tm <- max_tm
  res_melt_ACov_fit_LL.EP$step_size_td <- step_size_td
  res_melt_ACov_fit_LL.EP$step_size_tm <- step_size_tm
  res_melt_ACov_fit_LL.EP$JKchrom <- JKchrom
  
  
  
  
  
  ### plot all together
  EP_diff_values_all_data <- rbind(res_melt_ACov_fit_LL.EP)
  EP_fit_list[[x]] <- EP_diff_values_all_data
  
  
  
  
}
EP_fit_df <- do.call("rbind",EP_fit_list)

write.csv(x = Single_Sample_SP_ACov,file = as.character(output_SP),quote = F,row.names = F)
write.csv(x = Single_Sample_EP_ACov,file = as.character(output_EP),quote = F,row.names = F)
write.csv(x = EP_fit_df,file = as.character(output_EP_grid),quote = F,row.names = F)

Single_Sample_SP_summary <- Single_Sample_SP_ACov %>% mutate(td = NA )%>% 
  dplyr::select(tm_SP,td,sample,LL_SP) %>% dplyr::rename(tm = tm_SP,LL =  LL_SP) %>% mutate(method = "SP")

Single_Sample_EP_summary <- Single_Sample_EP_ACov %>% 
  dplyr::select(tm_EP,td_EP,sample,LL_EP) %>% dplyr::rename(tm = tm_EP,td = td_EP,LL =  LL_EP) %>% mutate(method = "EP_nls")

EP_grid_fit_summary <- EP_fit_df %>% group_by(ind) %>% filter(value0 == 0) %>%
  dplyr::select(Var2,Var1,ind,value) %>% dplyr::rename(tm = Var2, td = Var1, sample = ind, LL =  value) %>% mutate(method = "EP_grid")

JK_summary <- rbind(Single_Sample_SP_summary,Single_Sample_EP_summary,EP_grid_fit_summary)

write.csv(x = JK_summary,file = as.character(output_summary),quote = F,row.names = F)
