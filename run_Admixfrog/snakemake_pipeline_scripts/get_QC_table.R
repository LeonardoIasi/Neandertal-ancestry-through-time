library(ggplot2)
library(ggrepel)
library(tidyverse)
library(egg)


Meta_table <- readr::read_csv(snakemake@input$Meta_table)
Contamination_table <- read.table(snakemake@input$Contamination_table,header=T)
Posterior_Ancestry_estimates_table <- read.table(snakemake@input$Posterior_Ancestry_estimates_table,header = T)
Coverage_table = read.table(snakemake@input$Coverage_table,header=T)
f4_ratios_table <- read_csv(snakemake@input$f4_ratios_table)
EMH = as.logical(snakemake@params$EMH)


if(EMH == T){
  f4_ratios_table_output <- read_csv(snakemake@input$f4_ratios_table_output)
}

QC_table=inner_join(Contamination_table,Coverage_table, by = "sample_name") %>% 
  mutate(sample_name = gsub(pattern = "Old2200k",replacement = "", sample_name) ) %>%
  mutate(sample_name = gsub(pattern = "2200k",replacement = "",sample_name) ) %>%
  mutate(sample_name = gsub(pattern = "1240k",replacement = "",sample_name) )


if(EMH == T){
  Meta_table = Meta_table %>% select(c(admixfrog_sample_name,Country,World_Region,SGDP_Superpop,Latitude,Longitude,`Mean Calibrated Date OxCal 4.4.4  IntCal20 in BP`,
  `Calibrated Date OxCal 4.4.4  IntCal20 in BP 95.5 % lower`,`Calibrated Date OxCal 4.4.4  IntCal20 in BP 95.5 % higher`,
  `Most Likely Date Lower`,`Most Likely Date Higher`,`Most Likely Mean Date`,Average_DP,Data_source,Sex)) %>%
    mutate(pop = NA) %>% mutate(time = "ancient") 
  #Meta_table$admixfrog_sample_name[Meta_table$admixfrog_sample_name == "SalkhitArchAdm"] <- "Salkhit"
  QC_table=merge(QC_table,Meta_table, by.x = "sample_name", by.y = "admixfrog_sample_name")
  QC_table = QC_table %>%
    dplyr::rename(
      BP_mean = `Mean Calibrated Date OxCal 4.4.4  IntCal20 in BP`,
      BP_mean_lower = `Calibrated Date OxCal 4.4.4  IntCal20 in BP 95.5 % lower`,
      BP_mean_higher = `Calibrated Date OxCal 4.4.4  IntCal20 in BP 95.5 % higher`,
      ML_BP_Lower = `Most Likely Date Lower`,
      ML_BP_Higher = `Most Likely Date Higher`,
      ML_BP_Mean = `Most Likely Mean Date`,
      ave_cov_reported = Average_DP,
      ave_cov_array = ave_cov,
      world_region = World_Region
    ) %>% mutate(across(c(ave_cov_array), as.numeric)) %>% mutate(across(c(Data_source), as.factor))
} else {
  Meta_table = Meta_table %>% filter(dataset == "SGDP") %>%
    select(sample,pop,superpop,sex,Latitude,Longitude) %>% mutate(time = "present") %>% mutate(ave_cov = NA ) %>% mutate(Average_DP = 30) %>%
    mutate(Mean_date_BP = 0) %>% mutate(Data_source = "Genotypes")
  QC_table= QC_table %>% select(-ave_cov) %>% merge(.,Meta_table, by.x = "sample_name", by.y = "sample")
  QC_table = QC_table %>% mutate(ML_BP_Mean=Mean_date_BP) %>%
    dplyr::rename(
      BP_mean = Mean_date_BP,
      ave_cov_reported = Average_DP,
      ave_cov_array = ave_cov
    ) %>% mutate(across(c(BP_mean,ML_BP_Mean,ave_cov_reported,ave_cov_array), as.numeric)) %>% mutate(across(c(Data_source), as.factor))
}



# long table
f4_ratios_long <- f4_ratios_table %>% dplyr::rename(estimate = f4r,sample_name = sample) %>% mutate(ancestry_method="f4r") %>%
  mutate(lower=NA) %>% mutate(upper=NA) %>% 
  mutate(sample_name = gsub(pattern = "Old2200k",replacement = "", sample_name) ) %>%
  mutate(sample_name = gsub(pattern = "2200k",replacement = "",sample_name) ) %>%
  mutate(sample_name = gsub(pattern = "1240k",replacement = "",sample_name) )
if(EMH == T){
  f4_ratios_long_output <- f4_ratios_table_output %>% dplyr::rename(estimate = f4r,sample_name = sample) %>% mutate(ancestry_method="f4r_gtLH") %>%
    mutate(lower=NA) %>% mutate(upper=NA) 
  f4_ratios_long = rbind(f4_ratios_long,f4_ratios_long_output) %>% 
    mutate(sample_name = gsub(pattern = "Old2200k",replacement = "", sample_name) ) %>%
    mutate(sample_name = gsub(pattern = "2200k",replacement = "",sample_name) ) %>%
    mutate(sample_name = gsub(pattern = "1240k",replacement = "",sample_name) )
}
PP_ancestry_states <- Posterior_Ancestry_estimates_table %>% filter(state!="AFR") %>% 
  select(sample,state,mean,lower,upper) %>% dplyr::rename(sample_name=sample,Ancestry=state,estimate = mean) %>% mutate(data="all",ancestry_method="admixfrog") %>% 
  mutate(sample_name = gsub(pattern = "Old2200k",replacement = "", sample_name) ) %>%
  mutate(sample_name = gsub(pattern = "2200k",replacement = "",sample_name) ) %>%
  mutate(sample_name = gsub(pattern = "1240k",replacement = "",sample_name) )

QC_table_long <- rbind(f4_ratios_long,PP_ancestry_states) %>% inner_join(.,QC_table,by="sample_name")

QC_table_wide=QC_table_long %>%
  pivot_wider(
    names_from = c(data, ancestry_method,Ancestry),
    values_from = c(estimate,lower,upper)
  )


write.csv(QC_table_wide,snakemake@output$QC_table_wide,row.names = F)
write.csv(QC_table_long,snakemake@output$QC_table_long,row.names = F)
