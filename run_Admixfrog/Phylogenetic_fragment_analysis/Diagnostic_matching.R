library(tidyverse)
rle_path <- snakemake@input$input_anno_segment_file
diagnostic_pileup_path <- snakemake@input$diagnostic_pileup_file
diagnostic_site_ref_path <- snakemake@input$diagnostic_ref
bed_file <- snakemake@input$bed_file
DiagSites <- as.character(snakemake@wildcards$DiagSites)
type_ = as.character(snakemake@params$filter_type)
target_ = as.character(snakemake@params$filter_target)
min_len = as.numeric(snakemake@params$filter_map_len)
n_bootstrap = as.numeric(snakemake@params$n_bootstrap)
output_raw = snakemake@output$diagnostic_match_file_raw
output_summary = snakemake@output$diagnostic_match_file_sum
output_summary_TS_only = snakemake@output$diagnostic_match_file_sum_TS_only

chrom_ <- c("1" ,"2" , "3" , "4",  "5" , "6"  ,"7" , "8" , "9" , "10" ,"11", "12" ,"13","14" ,"15", "16" ,"17" ,"18", "19" , "20" ,"21", "22" , "X" )

called_fragments <- read_csv(rle_path,col_types = list(chrom=col_character())) %>%
  filter(type %in% type_)  %>%
  filter(target %in% target_) %>% filter(map_len >= min_len)

diagnpstic_pileup <- read.csv(diagnostic_pileup_path,sep="\t")
diagnpstic_pileup$chrom <- as.character(diagnpstic_pileup$chrom)

diagnostic_sites_ref <- read.table(diagnostic_site_ref_path,
                                   col.names = c("chrom","pos","ancestral","derived","Mbuti","Vindija33.19","Altai","Chagyrskaya","Denisova","derived_groups"))
diagnostic_sites_ref$chrom <- as.character(diagnostic_sites_ref$chrom)

if(DiagSites != "all"){
  diagnostic_sites_ref <-  diagnostic_sites_ref %>% filter(derived_groups != "Altai,Chagyrskaya,Denisova,Mbuti,Vindija33.19")
}

match_table <- inner_join(diagnostic_sites_ref,diagnpstic_pileup,by=c("chrom"="chrom","pos"="pos.end"))
match_table$total_n_reads <- rowSums(match_table[,c("A","C","G","T")],na.rm = T)
match_table <- match_table %>% filter(total_n_reads > 0)

if(DiagSites == "MatchingAscertainment"){
  print("subseting sites")
  bed_ <- read.table(bed_file)
  match_table <- inner_join(match_table,bed_,by=c("chrom"="V1","pos"="V3"))
}

if(length(called_fragments$chrom) == 0){
  
  match_table = match_table %>% mutate(proportion_match_Mbuti = NA,
                                   proportion_match_Vindija33.19 = NA,
                                   proportion_match_Altai = NA,
                                   proportion_match_Chagyrskaya = NA,
                                   proportion_match_Denisova = NA,
                                   p_derived = NA,
                                   sample = NA,
                                   Mbuti_matches = NA,
                                   Vindija33.19_matches = NA,
                                   Altai_matches = NA,
                                   Chagyrskaya_matches = NA,
                                   Denisova_matches = NA,
                                   Mbuti_matches_derived = NA,
                                   Vindija33.19_matches_derived = NA,
                                   Altai_matches_derived = NA,
                                   Chagyrskaya_matches_derived = NA,
                                   Denisova_matches_derived = NA)
  
  output_called_fragments = called_fragments %>% mutate(chrom = NA,
                                                        start = NA,
                                                        end = NA,
                                                        score = NA,
                                                        target = NA,
                                                        type = NA,
                                                        map = NA,
                                                        pos = NA,
                                                        id = NA,
                                                        map_end = NA,
                                                        pos_end = NA,
                                                        id_end = NA,
                                                        len = NA,
                                                        map_len = NA,
                                                        pos_len = NA,
                                                        nscore = NA,
                                                        n_all_snps = NA,
                                                        all_n_AFR = NA,
                                                        all_n_NEA = NA,
                                                        all_n_DEN = NA,
                                                        frag_ID = NA,
                                                        all_reads = NA,
                                                        prop_matching_Mbuti = NA,
                                                        prop_matching_Vindija33.19 = NA,
                                                        prop_matching_Altai = NA,
                                                        prop_matching_Chagyrskaya = NA,
                                                        prop_matching_Denisova = NA,
                                                        prop_pHcount_matching_Mbuti = NA,
                                                        count_matching_Mbuti = NA,
                                                        count_Mbuti = NA,
                                                        prop_pHcount_matching_Vindija33.19 = NA,
                                                        count_matching_Vindija33.19 = NA,
                                                        count_Vindija33.19 = NA,
                                                        prop_pHcount_matching_Altai = NA,
                                                        count_matching_Altai = NA,
                                                        count_Altai = NA,
                                                        prop_pHcount_matching_Chagyrskaya = NA,
                                                        count_matching_Chagyrskaya = NA,
                                                        count_Chagyrskaya = NA,
                                                        prop_pHcount_matching_Denisova = NA,
                                                        count_matching_Denisova = NA,
                                                        count_Denisova = NA,
                                                        BT_prop_pHcount_matching_Mbuti_mean = NA,
                                                        BT_prop_pHcount_matching_Mbuti_min = NA,
                                                        BT_prop_pHcount_matching_Mbuti_max = NA,
                                                        BT_prop_pHcount_matching_Mbuti_sd = NA,
                                                        BT_prop_pHcount_matching_Vindija33.19_mean = NA,
                                                        BT_prop_pHcount_matching_Vindija33.19_min = NA,
                                                        BT_prop_pHcount_matching_Vindija33.19_max = NA,
                                                        BT_prop_pHcount_matching_Vindija33.19_sd = NA,
                                                        BT_prop_pHcount_matching_Altai_mean = NA,
                                                        BT_prop_pHcount_matching_Altai_min = NA,
                                                        BT_prop_pHcount_matching_Altai_max = NA,
                                                        BT_prop_pHcount_matching_Altai_sd = NA,
                                                        BT_prop_pHcount_matching_Chagyrskaya_mean = NA,
                                                        BT_prop_pHcount_matching_Chagyrskaya_min = NA,
                                                        BT_prop_pHcount_matching_Chagyrskaya_max = NA,
                                                        BT_prop_pHcount_matching_Chagyrskaya_sd = NA,
                                                        BT_prop_pHcount_matching_Denisova_mean = NA,
                                                        BT_prop_pHcount_matching_Denisova_min = NA,
                                                        BT_prop_pHcount_matching_Denisova_max = NA,
                                                        BT_prop_pHcount_matching_Denisova_sd = NA)
  
     write.table(x=match_table,file = output_raw,quote = F,row.names = F,sep="\t")
     write.table(x=output_called_fragments,file = output_summary,quote = F,row.names = F,sep="\t")
     write.table(x=output_called_fragments,file = output_summary_TS_only,quote = F,row.names = F,sep="\t")
} else {


  assigne_derived <- function(data,read_cutoff) {
    ref_allele <- data["derived"]
    matching_reads <- as.numeric(data[paste0(ref_allele)])
    ifelse(matching_reads >= read_cutoff ,1,0)
  }
  
  assigne_derived_pseudoHaploid <- function(data) {
    ref_allele <- data["derived"]
    matching_reads <- as.numeric(data[paste0(ref_allele)])
    matching_reads/as.numeric(data["total_n_reads"])
  }
  
  # Add columns that evaluates the sample to be derived or ancestral
  match_table$p_derived <- apply(match_table, 1, assigne_derived_pseudoHaploid)
  match_table$p_derived[is.nan(match_table$p_derived)]<-NA
  match_table$sample <- rbinom(length(match_table$p_derived),1,match_table$p_derived)
  


  # count all matching reads
  count_matching_reads <- function(data, reference_column) {
    ref_allele <- ifelse(data[[reference_column]] == 1, data["derived"], data["ancestral"])
    matching_reads <- as.numeric(data[paste0(ref_allele)])
    sum(matching_reads, na.rm = TRUE)
  }
  
  # count all matching reads
  count_matching_derived_reads <- function(data, reference_column) {
    if(data[[reference_column]] == 0){
      return(0)
    } else{
      ref_allele <- data["derived"]
      matching_reads <- as.numeric(data[paste0(ref_allele)])
      return(sum(matching_reads, na.rm = TRUE))
    }

  }

  reference_columns = c("Mbuti","Vindija33.19","Altai","Chagyrskaya","Denisova")

  # Add columns for the count of matching reads for each reference individual
  for (ref_col in reference_columns) {
    match_table[[paste0(ref_col, "_matches")]] <- apply(match_table, 1, count_matching_reads, reference_column = ref_col)
  }
  
  for (ref_col in reference_columns) {
    match_table[[paste0(ref_col, "_matches_derived")]] <- apply(match_table, 1, count_matching_derived_reads, reference_column = ref_col)
  }

  
  #prop_read_per_seg = match_table %>% filter(total_n_reads >= 1) %>% group_by(ID) %>% summarize(all_reads = sum(total_n_reads),
  #                                                                                              prop_matching_Mbuti = (sum(Mbuti_matches_derived,na.rm = T)/sum(total_n_reads)),
  #                                                                                              prop_matching_Vindija33.19 = (sum(Vindija33.19_matches_derived,na.rm = T)/sum(total_n_reads)),
  #                                                                                              prop_matching_Altai = (sum(Altai_matches_derived,na.rm = T)/sum(total_n_reads)),
  #                                                                                              prop_matching_Chagyrskaya = (sum(Chagyrskaya_matches_derived,na.rm = T)/sum(total_n_reads)),
  #                                                                                              prop_matching_Denisova = (sum(Denisova_matches_derived,na.rm = T)/sum(total_n_reads)))
  
  prop_read_per_seg = match_table %>% mutate(n_derived_reads = total_n_reads * p_derived) %>% filter(n_derived_reads >= 1) %>% group_by(ID) %>% summarize(all_reads = sum(n_derived_reads),
                                                                                                prop_matching_Mbuti = (sum(Mbuti_matches_derived,na.rm = T)/sum(n_derived_reads)),
                                                                                                prop_matching_Vindija33.19 = (sum(Vindija33.19_matches_derived,na.rm = T)/sum(n_derived_reads)),
                                                                                                prop_matching_Altai = (sum(Altai_matches_derived,na.rm = T)/sum(n_derived_reads)),
                                                                                                prop_matching_Chagyrskaya = (sum(Chagyrskaya_matches_derived,na.rm = T)/sum(n_derived_reads)),
                                                                                                prop_matching_Denisova = (sum(Denisova_matches_derived,na.rm = T)/sum(n_derived_reads)))
  

  prop_read_per_seg_TransversionOnly = match_table %>% 
    mutate(n_derived_reads = total_n_reads * p_derived) %>% filter(n_derived_reads >= 1) %>% 
    filter(!(ancestral=="C" & derived == "T")  & !(ancestral=="T" & derived == "C") & !(ancestral=="A" & derived == "G") & !(ancestral=="G" & derived == "A") ) %>%
    group_by(ID) %>% 
    summarize(all_reads = sum(n_derived_reads),
              prop_matching_Mbuti = (sum(Mbuti_matches_derived,na.rm = T)/sum(n_derived_reads)),
              prop_matching_Vindija33.19 = (sum(Vindija33.19_matches_derived,na.rm = T)/sum(n_derived_reads)),
              prop_matching_Altai = (sum(Altai_matches_derived,na.rm = T)/sum(n_derived_reads)),
              prop_matching_Chagyrskaya = (sum(Chagyrskaya_matches_derived,na.rm = T)/sum(n_derived_reads)),
              prop_matching_Denisova = (sum(Denisova_matches_derived,na.rm = T)/sum(n_derived_reads)))
  
                                                                                                     
  #count_pseudoHapl_derived_per_seg <- match_table%>% filter(total_n_reads >= 1) %>% group_by(ID) %>% summarize(
  #  prop_pHcount_matching_Mbuti = (sample %*% Mbuti) / n(),
  #  count_matching_Mbuti = (sample %*% Mbuti),
  #  count_Mbuti = sum(Mbuti),
  #  prop_pHcount_matching_Vindija33.19 = (sample %*% Vindija33.19)  / n(),
  #  count_matching_Vindija33.19 = (sample %*% Vindija33.19),
  #  count_Vindija33.19 = sum(Vindija33.19),
  #  prop_pHcount_matching_Altai = (sample %*% Altai)  / n(),
  #  count_matching_Altai = (sample %*% Altai),
  #  count_Altai = sum(Altai),
  #  prop_pHcount_matching_Chagyrskaya = (sample %*% Chagyrskaya)  / n(),
  #  count_matching_Chagyrskaya = (sample %*% Chagyrskaya),
  #  count_Chagyrskaya = sum(Chagyrskaya),
  #  prop_pHcount_matching_Denisova = (sample %*% Denisova)  / n(),
  #  count_matching_Denisova = (sample %*% Denisova),
  #  count_Denisova = sum(Denisova))
  
  count_pseudoHapl_derived_per_seg <- match_table%>% filter(total_n_reads >= 1) %>% group_by(ID) %>% summarize(
    prop_pHcount_matching_Mbuti = (sample %*% Mbuti) / sum(sample),
    count_matching_Mbuti = (sample %*% Mbuti),
    count_Mbuti = sum(Mbuti),
    prop_pHcount_matching_Vindija33.19 = (sample %*% Vindija33.19)  / sum(sample),
    count_matching_Vindija33.19 = (sample %*% Vindija33.19),
    count_Vindija33.19 = sum(Vindija33.19),
    prop_pHcount_matching_Altai = (sample %*% Altai)  / sum(sample),
    count_matching_Altai = (sample %*% Altai),
    count_Altai = sum(Altai),
    prop_pHcount_matching_Chagyrskaya = (sample %*% Chagyrskaya)  / sum(sample),
    count_matching_Chagyrskaya = (sample %*% Chagyrskaya),
    count_Chagyrskaya = sum(Chagyrskaya),
    prop_pHcount_matching_Denisova = (sample %*% Denisova)  / sum(sample),
    count_matching_Denisova = (sample %*% Denisova),
    count_Denisova = sum(Denisova))
  
  count_pseudoHapl_derived_per_seg_TransversionOnly <- match_table %>% filter(total_n_reads >= 1) %>% 
    filter(!(ancestral=="C" & derived == "T")  & !(ancestral=="T" & derived == "C") & !(ancestral=="A" & derived == "G") & !(ancestral=="G" & derived == "A") ) %>%
    group_by(ID) %>% summarize(
      prop_pHcount_matching_Mbuti = (sample %*% Mbuti) / sum(sample),
      count_matching_Mbuti = (sample %*% Mbuti),
      count_Mbuti = sum(Mbuti),
      prop_pHcount_matching_Vindija33.19 = (sample %*% Vindija33.19)  / sum(sample),
      count_matching_Vindija33.19 = (sample %*% Vindija33.19),
      count_Vindija33.19 = sum(Vindija33.19),
      prop_pHcount_matching_Altai = (sample %*% Altai)  / sum(sample),
      count_matching_Altai = (sample %*% Altai),
      count_Altai = sum(Altai),
      prop_pHcount_matching_Chagyrskaya = (sample %*% Chagyrskaya)  / sum(sample),
      count_matching_Chagyrskaya = (sample %*% Chagyrskaya),
      count_Chagyrskaya = sum(Chagyrskaya),
      prop_pHcount_matching_Denisova = (sample %*% Denisova)  / sum(sample),
      count_matching_Denisova = (sample %*% Denisova),
      count_Denisova = sum(Denisova))
  
  
  Bootstrap_prop_pHcount_matching_Mbuti <- matrix(NA,ncol = n_bootstrap,nrow = length(unique(count_pseudoHapl_derived_per_seg$ID)))
  Bootstrap_prop_pHcount_matching_Vindija33.19 <- matrix(NA,ncol = n_bootstrap,nrow = length(unique(count_pseudoHapl_derived_per_seg$ID)))
  Bootstrap_prop_pHcount_matching_Altai <- matrix(NA,ncol = n_bootstrap,nrow = length(unique(count_pseudoHapl_derived_per_seg$ID)))
  Bootstrap_prop_pHcount_matching_Chagyrskaya <- matrix(NA,ncol = n_bootstrap,nrow = length(unique(count_pseudoHapl_derived_per_seg$ID)))
  Bootstrap_prop_pHcount_matching_Denisova <- matrix(NA,ncol = n_bootstrap,nrow = length(unique(count_pseudoHapl_derived_per_seg$ID)))
  
  Bootstrap_prop_pHcount_matching_Mbuti_TS <- matrix(NA,ncol = n_bootstrap,nrow = length(unique(count_pseudoHapl_derived_per_seg_TransversionOnly$ID)))
  Bootstrap_prop_pHcount_matching_Vindija33.19_TS <- matrix(NA,ncol = n_bootstrap,nrow = length(unique(count_pseudoHapl_derived_per_seg_TransversionOnly$ID)))
  Bootstrap_prop_pHcount_matching_Altai_TS <- matrix(NA,ncol = n_bootstrap,nrow = length(unique(count_pseudoHapl_derived_per_seg_TransversionOnly$ID)))
  Bootstrap_prop_pHcount_matching_Chagyrskaya_TS <- matrix(NA,ncol = n_bootstrap,nrow = length(unique(count_pseudoHapl_derived_per_seg_TransversionOnly$ID)))
  Bootstrap_prop_pHcount_matching_Denisova_TS <- matrix(NA,ncol = n_bootstrap,nrow = length(unique(count_pseudoHapl_derived_per_seg_TransversionOnly$ID)))
  
  for(i in 1:n_bootstrap){
    
    bootstrap_x <- match_table %>% filter(total_n_reads >= 1) %>% 
      mutate(sample_bootstrap = rbinom(length(p_derived),1,p_derived)) %>% 
      group_by(ID) %>% summarize(
      prop_pHcount_matching_Mbuti = (sample_bootstrap %*% Mbuti) / sum(sample_bootstrap),
      prop_pHcount_matching_Vindija33.19 = (sample_bootstrap %*% Vindija33.19)  / sum(sample_bootstrap),
      prop_pHcount_matching_Altai = (sample_bootstrap %*% Altai)  / sum(sample_bootstrap),
      prop_pHcount_matching_Chagyrskaya = (sample_bootstrap %*% Chagyrskaya)  / sum(sample_bootstrap),
      prop_pHcount_matching_Denisova = (sample_bootstrap %*% Denisova)  / sum(sample_bootstrap))
    
    Bootstrap_prop_pHcount_matching_Mbuti[,i] <- bootstrap_x$prop_pHcount_matching_Mbuti
    Bootstrap_prop_pHcount_matching_Vindija33.19[,i] <- bootstrap_x$prop_pHcount_matching_Vindija33.19
    Bootstrap_prop_pHcount_matching_Altai[,i] <- bootstrap_x$prop_pHcount_matching_Altai
    Bootstrap_prop_pHcount_matching_Chagyrskaya[,i] <- bootstrap_x$prop_pHcount_matching_Chagyrskaya
    Bootstrap_prop_pHcount_matching_Denisova[,i] <- bootstrap_x$prop_pHcount_matching_Denisova
    
    
    bootstrap_TS_x <- match_table %>% filter(total_n_reads >= 1) %>% 
      filter(!(ancestral=="C" & derived == "T")  & !(ancestral=="T" & derived == "C") & !(ancestral=="A" & derived == "G") & !(ancestral=="G" & derived == "A") ) %>%
      mutate(sample_bootstrap_TS = rbinom(length(p_derived),1,p_derived)) %>% 
      group_by(ID) %>% summarize(
        prop_pHcount_matching_Mbuti = (sample_bootstrap_TS %*% Mbuti) / sum(sample_bootstrap_TS),
        prop_pHcount_matching_Vindija33.19 = (sample_bootstrap_TS %*% Vindija33.19)  / sum(sample_bootstrap_TS),
        prop_pHcount_matching_Altai = (sample_bootstrap_TS %*% Altai)  / sum(sample_bootstrap_TS),
        prop_pHcount_matching_Chagyrskaya = (sample_bootstrap_TS %*% Chagyrskaya)  / sum(sample_bootstrap_TS),
        prop_pHcount_matching_Denisova = (sample_bootstrap_TS %*% Denisova) / sum(sample_bootstrap_TS))
    
    Bootstrap_prop_pHcount_matching_Mbuti_TS[,i] <- bootstrap_TS_x$prop_pHcount_matching_Mbuti
    Bootstrap_prop_pHcount_matching_Vindija33.19_TS[,i] <- bootstrap_TS_x$prop_pHcount_matching_Vindija33.19
    Bootstrap_prop_pHcount_matching_Altai_TS[,i] <- bootstrap_TS_x$prop_pHcount_matching_Altai
    Bootstrap_prop_pHcount_matching_Chagyrskaya_TS[,i] <- bootstrap_TS_x$prop_pHcount_matching_Chagyrskaya
    Bootstrap_prop_pHcount_matching_Denisova_TS[,i] <- bootstrap_TS_x$prop_pHcount_matching_Denisova
    
    
  }

  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Mbuti_mean <- apply(Bootstrap_prop_pHcount_matching_Mbuti, 1, FUN = mean)  
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Mbuti_min <- apply(Bootstrap_prop_pHcount_matching_Mbuti, 1, FUN = min)
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Mbuti_max <- apply(Bootstrap_prop_pHcount_matching_Mbuti, 1, FUN = max)
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Mbuti_sd <- apply(Bootstrap_prop_pHcount_matching_Mbuti, 1, FUN = sd)

  
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Vindija33.19_mean <- apply(Bootstrap_prop_pHcount_matching_Vindija33.19, 1, FUN = mean)   
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Vindija33.19_min <- apply(Bootstrap_prop_pHcount_matching_Vindija33.19, 1, FUN = min)
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Vindija33.19_max <- apply(Bootstrap_prop_pHcount_matching_Vindija33.19, 1, FUN = max)
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Vindija33.19_sd <- apply(Bootstrap_prop_pHcount_matching_Vindija33.19, 1, FUN = sd)
  
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Altai_mean <- apply(Bootstrap_prop_pHcount_matching_Altai, 1, FUN = mean)  
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Altai_min <- apply(Bootstrap_prop_pHcount_matching_Altai, 1, FUN = min)
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Altai_max <- apply(Bootstrap_prop_pHcount_matching_Altai, 1, FUN = max)
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Altai_sd <- apply(Bootstrap_prop_pHcount_matching_Altai, 1, FUN = sd)
  
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Chagyrskaya_mean <- apply(Bootstrap_prop_pHcount_matching_Chagyrskaya, 1, FUN = mean)
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Chagyrskaya_min <- apply(Bootstrap_prop_pHcount_matching_Chagyrskaya, 1, FUN = min)
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Chagyrskaya_max <- apply(Bootstrap_prop_pHcount_matching_Chagyrskaya, 1, FUN = max)
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Chagyrskaya_sd <- apply(Bootstrap_prop_pHcount_matching_Chagyrskaya, 1, FUN = sd)
  
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Denisova_mean <- apply(Bootstrap_prop_pHcount_matching_Denisova, 1, FUN = mean)
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Denisova_min <- apply(Bootstrap_prop_pHcount_matching_Denisova, 1, FUN = min)
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Denisova_max <- apply(Bootstrap_prop_pHcount_matching_Denisova, 1, FUN = max)
  count_pseudoHapl_derived_per_seg$BT_prop_pHcount_matching_Denisova_sd <- apply(Bootstrap_prop_pHcount_matching_Denisova, 1, FUN = sd)
  
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Mbuti_mean <- apply(Bootstrap_prop_pHcount_matching_Mbuti_TS, 1, FUN = mean)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Mbuti_min <- apply(Bootstrap_prop_pHcount_matching_Mbuti_TS, 1, FUN = min)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Mbuti_max <- apply(Bootstrap_prop_pHcount_matching_Mbuti_TS, 1, FUN = max)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Mbuti_sd <- apply(Bootstrap_prop_pHcount_matching_Mbuti_TS, 1, FUN = sd)
  
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Vindija33.19_mean <- apply(Bootstrap_prop_pHcount_matching_Vindija33.19_TS, 1, FUN = mean)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Vindija33.19_min <- apply(Bootstrap_prop_pHcount_matching_Vindija33.19_TS, 1, FUN = min)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Vindija33.19_max <- apply(Bootstrap_prop_pHcount_matching_Vindija33.19_TS, 1, FUN = max)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Vindija33.19_sd <- apply(Bootstrap_prop_pHcount_matching_Vindija33.19_TS, 1, FUN = sd)
  
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Altai_mean <- apply(Bootstrap_prop_pHcount_matching_Altai_TS, 1, FUN = mean)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Altai_min <- apply(Bootstrap_prop_pHcount_matching_Altai_TS, 1, FUN = min)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Altai_max <- apply(Bootstrap_prop_pHcount_matching_Altai_TS, 1, FUN = max)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Altai_sd <- apply(Bootstrap_prop_pHcount_matching_Altai_TS, 1, FUN = sd)
  
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Chagyrskaya_mean <- apply(Bootstrap_prop_pHcount_matching_Chagyrskaya_TS, 1, FUN = mean)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Chagyrskaya_min <- apply(Bootstrap_prop_pHcount_matching_Chagyrskaya_TS, 1, FUN = min)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Chagyrskaya_max <- apply(Bootstrap_prop_pHcount_matching_Chagyrskaya_TS, 1, FUN = max)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Chagyrskaya_sd <- apply(Bootstrap_prop_pHcount_matching_Chagyrskaya_TS, 1, FUN = sd)
  
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Denisova_mean <- apply(Bootstrap_prop_pHcount_matching_Denisova_TS, 1, FUN = mean)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Denisova_min <- apply(Bootstrap_prop_pHcount_matching_Denisova_TS, 1, FUN = min)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Denisova_max <- apply(Bootstrap_prop_pHcount_matching_Denisova_TS, 1, FUN = max)
  count_pseudoHapl_derived_per_seg_TransversionOnly$BT_prop_pHcount_matching_Denisova_sd <- apply(Bootstrap_prop_pHcount_matching_Denisova_TS, 1, FUN = sd)

  match_table_summary <- inner_join(prop_read_per_seg,count_pseudoHapl_derived_per_seg,by="ID")
  match_table_summary_transversions_only <- inner_join(prop_read_per_seg_TransversionOnly,count_pseudoHapl_derived_per_seg_TransversionOnly,by="ID")

  output_called_fragments <- left_join(called_fragments,match_table_summary,by=c("frag_ID"="ID")) 
  output_called_fragments_transversions_only <- left_join(called_fragments,match_table_summary_transversions_only,by=c("frag_ID"="ID")) 

  write.table(x=match_table,file = output_raw,quote = F,row.names = F,sep="\t")
  write.table(x=output_called_fragments,file = output_summary,quote = F,row.names = F,sep="\t")
  write.table(x=output_called_fragments_transversions_only,file = output_summary_TS_only,quote = F,row.names = F,sep="\t")

}
