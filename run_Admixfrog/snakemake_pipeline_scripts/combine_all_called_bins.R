library(tidyverse)
bin_file_path_list = snakemake@input$res_bin_call
ref_bin_file_path = snakemake@input$bin_ref_file
min_len = snakemake@wildcards$bin_called_min_len
min_len_pos = snakemake@wildcards$bin_called_min_len_pos
min_n_all_SNP = snakemake@wildcards$allSNP
min_n_SNP = snakemake@wildcards$SNP
split_char = snakemake@params$split_char
Ancestry = snakemake@wildcards$called_Ancestry
All_Ancestries = snakemake@params$All_ancestries
output = snakemake@output$all_called_bins
output2 = snakemake@output$all_called_bins_frag_ID

ref_bin_file = read_csv(ref_bin_file_path,col_types = list(chrom=col_character())) %>% dplyr::select(chrom,map,pos,id)
ref_bin_file$map <- round(ref_bin_file$map,3)
ref_bins_frag_ID = ref_bin_file

ref_bins_fn <- function(bin_file_path,ref_bin,Ancestry,All_Ancestries,min_len,min_len_pos,split_char){
  sample_name = str_split(bin_file_path,"/")[[1]] %>% tail(.,n=1) %>% str_split(.,split_char)
  sample_name = sample_name[[1]][1]
  print(paste("processing ",sample_name,sep=""))
  bin_file = read.csv(bin_file_path)
  bin_file$chrom <- as.character(bin_file$chrom)
  bin_file$map <- round(bin_file$map,3)
  ind_bins = bin_file %>% filter(called == Ancestry) %>% filter(total_map_len>= min_len,total_pos_len >= min_len_pos,n_all_snps >= min_n_all_SNP,n_snps >= min_n_SNP) 
  if(length(ind_bins$chrom) == 0){
    ref_bin_file[,sample_name] <- 0 
  } else {
    ind_bins$bin <- rowSums(ind_bins[All_Ancestries[str_detect(All_Ancestries,Ancestry)]])
    
    xx = ind_bins %>% dplyr::select(chrom,map,bin) %>% full_join(ref_bin_file, ., by = c("chrom","map"))
    xx$bin[is.na(xx$bin)] = 0
    ref_bin_file[,sample_name] <- xx %>% dplyr::select(chrom,map,bin) %>% inner_join(.,ref_bin_file) %>% distinct() %>% dplyr::select(bin)
  }

  return(ref_bin_file)
}

for(i in bin_file_path_list){
  ref_bin_file = ref_bins_fn(bin_file_path = i,ref_bin=ref_bin_file, Ancestry=Ancestry,All_Ancestries=All_Ancestries,min_len=min_len,min_len_pos=min_len_pos,split_char=split_char)
}

write.csv(file = output,x = ref_bin_file,quote = F,row.names = F)

ref_bins_frag_ID_fn <- function(bin_file_path,ref_bin,Ancestry,All_Ancestries,min_len,min_len_pos,split_char){
  sample_name = str_split(bin_file_path,"/")[[1]] %>% tail(.,n=1) %>% str_split(.,split_char)
  sample_name = sample_name[[1]][1]
  print(paste("processing ",sample_name,sep=""))
  bin_file = read.csv(bin_file_path)
  bin_file$chrom <- as.character(bin_file$chrom)
  bin_file$map <- round(bin_file$map,3)
  ind_bins = bin_file %>% filter(called == Ancestry) %>% filter(total_map_len>= min_len,total_pos_len >= min_len_pos,n_all_snps >= min_n_all_SNP,n_snps >= min_n_SNP) 
  if(length(ind_bins$chrom) == 0){
    ref_bin[,sample_name] <- NA 
  } else {
    xx = ind_bins %>% dplyr::select(chrom,map,frag_ID) %>% full_join(ref_bin, ., by = c("chrom","map"))
    ref_bin[,sample_name] <- xx %>% dplyr::select(chrom,map,frag_ID) %>% inner_join(.,ref_bin) %>% distinct() %>% dplyr::select(frag_ID)
  }
  
  return(ref_bin)
}

for(i in bin_file_path_list){
  ref_bins_frag_ID = ref_bins_frag_ID_fn(i,ref_bin=ref_bins_frag_ID, Ancestry=Ancestry,All_Ancestries=All_Ancestries,min_len=min_len,min_len_pos=min_len_pos,split_char=split_char)
}

write.csv(file = output2,x = ref_bins_frag_ID,quote = F,row.names = F)
