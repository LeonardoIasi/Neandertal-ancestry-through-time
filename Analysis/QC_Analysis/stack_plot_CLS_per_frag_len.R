library(tidyverse)
obs_path_list = snakemake@input$p_comp_table_data
meta_data = snakemake@input$meta_data
ind_ = snakemake@wildcards$sample
TP_ = snakemake@wildcards$TruePenalty
OP_ = snakemake@wildcards$penalty
minLen = snakemake@params$minLen
maxLen = snakemake@params$maxLen
plot_var = snakemake@params$plot_var
plot_asc = snakemake@params$plot_asc
output = snakemake@output$p_frag_stats
output2 = snakemake@output$p_frag_stats_RData

print(output)
print(output2)

res_obs = do.call(rbind, lapply(obs_path_list, read.table, header=TRUE, sep=","))

meta_data = readr::read_csv(meta_data)

ave_cov = meta_data %>% filter(sample_name == ind_) %>% dplyr::select(ave_cov_array)

res_obs %>%  filter(Test_status != "true") %>% mutate(map_len2=ifelse(!is.na(map_len),map_len,est_map_len))  %>% filter(map_len2 <= maxLen)  %>%
  filter(TruePenalty == TP_) %>% filter(Obs_penalty == OP_) %>% 
  mutate(map_len_bin = cut(map_len2, breaks = 30),
         map_len_bin_temp = gsub(pattern = "\\(|\\[|\\)|\\]", replacement = "", map_len_bin)) %>% 
  separate(col = map_len_bin_temp,into= c("lwr","upr"),sep=",") %>%
  mutate(map_len_bin_mean = (as.numeric(lwr)+as.numeric(upr))/2)-> Plot_data


#if(!is.null(output)){
#  if(plot_var == "AdmixFParams"){
#    print(paste("filtering for ",plot_asc," ascertainment"))
#    Plot_data = Plot_data %>% filter(Obs_Ascertainment == plot_asc)
#  }
#  q
#  P2 <- ggplot(Plot_data, aes(fill=CLS, y=map_len2,x=map_len_bin_mean)) + 
#    geom_bar(position="fill", stat="identity") +
#    facet_grid(vars(!! rlang::sym(plot_var) ), vars(as.factor(round(Obs_downsampling*ave_cov$ave_cov_array,2))),scales = "free")+
#    xlab("genetic length")+
#    ylab("proportion")+
#    ggtitle(paste(ind_,"True penalty: ",TP_," ; Observed penalty: ",OP_,sep = "")) 
#  
#  ggsave(filename = output,plot =  P2,device = "png")
#}


save.image(output2)



