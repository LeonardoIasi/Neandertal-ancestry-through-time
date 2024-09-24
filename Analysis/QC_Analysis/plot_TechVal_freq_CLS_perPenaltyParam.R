library(tidyverse)

obs_path_list = snakemake@input$p_comp_table_data

outfile = snakemake@output$p_frag_stats
outfile2 = paste(substr(outfile,1,nchar(outfile)-4),".csv",sep="")

res_obs = do.call(rbind, lapply(obs_path_list, read.table, header=TRUE, sep=","))

min_len = unique(res_obs$len_filter_all)
min_len_obs = unique(res_obs$len_filter_obs)
TP_ = unique(res_obs$TruePenalty)

write_csv(res_obs,outfile2)

res_obs %>%  filter(Test_status != "true") %>% group_by(Obs_penalty,Obs_Ascertainment,Obs_downsampling,CLS,AdmixFParams) %>% summarise(n=n()) %>%
  mutate(n = n) %>% with_groups(c(Obs_penalty,Obs_Ascertainment,Obs_downsampling,AdmixFParams),mutate,N=sum(n)) %>% mutate(freq=n/N) -> plot_data
plot_data$penalty <- as.factor(plot_data$Obs_penalty)

P1 <-   ggplot(plot_data, aes(x=Obs_Ascertainment,y= freq,color=penalty)) + 
  geom_point() +
  theme(axis.text.x=element_text(angle = -90)) +
  facet_grid(vars(CLS), vars(as.factor(Obs_downsampling)),scales = 'free')+
  xlab("") +
  ylab("") +
  ggtitle(paste("min len all/obs = ",min_len,"/",min_len+min_len_obs," cM; True penalty = ",TP_,sep="")) 

ggsave(outfile, P1,device = "png")