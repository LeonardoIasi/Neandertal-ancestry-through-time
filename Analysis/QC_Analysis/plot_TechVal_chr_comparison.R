library(tidyverse)

true_path = snakemake@input$true_comp
obs_path_list = snakemake@input$obs_comp
chr_true_path = snakemake@input$true_rle


min_len = snakemake@params$min_len
min_len_obs = snakemake@params$min_len_obs
asc_obs = snakemake@wildcards$snpset
type_ = snakemake@params$type
state_ = snakemake@params$state
chr_ = snakemake@params$chr
TP_ = snakemake@wildcards$TruePenalty
AF_ = snakemake@wildcards$admf_params


outfile2 = snakemake@output$p_comp_table_data
outfile3 = snakemake@output$p_comp_chr


res_obs_only = do.call(rbind, lapply(obs_path_list, read.table, header=TRUE, sep=","))
res_obs_only$Test_status <- "observation"
true = read.table(true_path,header = T,sep=",")
true$Test_status <- "true"

res_obs <- rbind(true,res_obs_only)
res_obs <- res_obs %>% select(., -est_target)
res_obs$abs_error_len <- abs(res_obs$err_start)+abs(res_obs$err_end)
res_obs$error_len <- res_obs$err_end - res_obs$err_start
res_obs$est_len <- res_obs$len + res_obs$error_len
res_obs$map_len <- res_obs$map_end - res_obs$map_start
res_obs$est_map_len <- res_obs$est_map_end - res_obs$est_map_start
res_obs$test_sample <- paste(res_obs$Obs_downsampling,res_obs$Obs_Ascertainment,sep="_")
res_obs$TruePenalty <- as.numeric(res_obs$True_penalty)
res_obs$TrueAscertainment <- as.numeric(res_obs$True_Ascertainment)
res_obs$AdmixFParams <- AF_

true_name <- true$True_sample[1]
print(true_name)

print(colnames(res_obs))

res_obs_filtered <- res_obs  %>% filter(Obs_sample != true_name & CLS=="TP" & map_len >= min_len_obs & est_map_len >= min_len_obs |Obs_sample != true_name & CLS=="FP" & est_map_len >= min_len_obs |Obs_sample != true_name & CLS=="GAP" & map_len >= min_len_obs & est_map_len >= min_len_obs |Obs_sample != true_name & CLS=="MERGED" & map_len >= min_len_obs & est_map_len >= min_len_obs |Obs_sample != true_name & CLS=="FN" & map_len >= min_len_obs | Obs_sample == true_name )
res_obs_filtered$len_filter_obs <- min_len_obs
res_obs_filtered$len_filter_all <- min_len
obs_penalty <- res_obs_filtered$Obs_penalty[res_obs_filtered$Test_status != "true"][1]

write_csv(res_obs_filtered,outfile2)


if(!is.null(outfile3)){
  source(paste(getwd(),"/plotting/lib_Leo.R",sep=""), chdir=T)
  chr_true <- read.csv(chr_true_path)
  chr_true = chr_true %>% filter(type == type_) 
  
  res_obs_filtered %>%  filter(chrom==chr_) %>% 
    ggplot() +
    bg_chrom(chr_true) +
    geom_tile(aes(x =(map_start +map_end)/ 2, width = map_end - map_start,
                  y=0.5, height=1,fill=CLS)) +
    facet_wrap(~test_sample, ncol=1, strip='left') +
    THEME +
    scale_alpha_continuous(range=c(0.3,1), limits=c(0.1, 2), name='Length(cM)') +
    scale_x_continuous(expand=c(0,0), name='Position (cM)') +
    scale_y_discrete(expand=c(0,0), name='State') +
    theme(strip.text.y = element_text(angle=180))+
    ggtitle(paste(asc_obs," penalty ",obs_penalty," chr ",chr_," min len all/obs = ",min_len,"/",min_len+min_len_obs," cM; True penalty ",TP_,sep="")) -> P2
  
  ggsave(outfile3, P2, width=20, height=11.5,device = "png")
}
