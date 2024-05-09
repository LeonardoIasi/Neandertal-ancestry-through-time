library(tidyverse)
bin_file_path = snakemake@input$res_bin
rle_file_path = snakemake@input$rle
ref_file_path = snakemake@input$ref
snp_file_path = snakemake@input$res_snp
target_ = snakemake@params$major_ref
gtmode = snakemake@params$gtmode
all = snakemake@params$type_alle
sample = snakemake@wildcards$sample
pops = strsplit(snakemake@wildcards$states,"-")[[1]]
output_bin = snakemake@output$res_bin_call
output_rle = snakemake@output$res_rle

# functions

'%notin%' <- Negate('%in%')

get_call_status <- function(rle_file,bin_file){
  bin_archaic = bin_file %>% filter(id %in% c(rle_file$start:rle_file$end )) %>% mutate(called = rle_file$target,
                                                                                        total_map_len = rle_file$map_len,
                                                                                        total_pos_len = rle_file$pos_len,
                                                                                        map_start = rle_file$map,
                                                                                        map_end = rle_file$map_end,
                                                                                        total_pos_len = rle_file$pos_len,
                                                                                        pos_start = rle_file$pos,
                                                                                        pos_end = rle_file$pos_end,
                                                                                        frag_ID = paste(sample,rle_file$target,rle_file$chrom,rle_file$start,rle_file$end,sep="_"))
  return(bin_archaic)
}

get_snp_status <- function(rle_file_x,bin_called_file,pops){
  rle_file_x = as.data.frame(rle_file_x)
  snp_rle_x = bin_called_file %>% filter(id %in% c(rle_file_x$start:rle_file_x$end )) %>% dplyr::summarise(n_snp = sum(n_snps,na.rm = T))
  rle_file_x["n_snp"] = snp_rle_x$n_snp
  for(pop in pops){
    varname_ = sprintf("n_%s", pop)
    new_varname_ = sprintf("all_n_%s", pop)
    snp_rle_x = bin_called_file %>% filter(id %in% c(rle_file_x$start:rle_file_x$end )) %>% dplyr::summarise(!!new_varname_ := sum(!!sym(varname_),na.rm = T))
    rle_file_x[new_varname_] = snp_rle_x[,new_varname_]
  }
  return(rle_file_x)
}

# read in data

bin_file = read_csv(bin_file_path,col_types = list(chrom=col_character()))

snp_file = read_csv(snp_file_path,col_types = list(chrom=col_character()))

if(gtmode == T){
  snp_file = snp_file %>% mutate(random_read = rbinom(length(snp_file$talt),1,snp_file$talt/(snp_file$talt + snp_file$tref)))
}

ref_file = read_csv(ref_file_path,col_types = list(chrom=col_character()))

if(any((paste0(pops,"_alt") %notin% colnames(ref_file)))){
  pop_missing = pops[(paste0(pops,"_alt") %notin% colnames(ref_file))]
  pops = pops[pops %notin% pop_missing]
  for(pop in pop_missing){
    pop_name = strsplit(pop,split = "=")[[1]][1]
    pop_inds = strsplit(strsplit(pop,split = "=")[[1]][2],split = "\\+")[[1]]
    ref_file[paste0(pop_name,"_alt")] <- rowSums(ref_file[,paste0(pop_inds,"_alt")])
    ref_file[paste0(pop_name,"_ref")] <- rowSums(ref_file[,paste0(pop_inds,"_ref")])
    pops <- c(pops,pop_name)
    
  }
}

if(all == T){
  print("Writing all types")
  rle_file = read_csv(rle_file_path,col_types = list(chrom=col_character()))  
  
}else{
  rle_file = read_csv(rle_file_path,col_types = list(chrom=col_character())) %>%
    filter(type == "state") 
}




bin_archaic = rle_file %>% filter(target != target_) %>% rowwise %>%
  do(X=get_call_status(rle_file=., bin_file=bin_file)) %>%
  unnest(X) %>% bind_cols(.)

bin_file$called = target_
bin_file$total_map_len = NA
bin_file$total_pos_len = NA
bin_file$map_start = NA
bin_file$map_end = NA
bin_file$pos_start = NA
bin_file$pos_end = NA
bin_file$frag_ID = NA


bin_called_file = rbind(bin_file[!bin_file$id %in% bin_archaic$id,],bin_archaic)

run_snp = snp_file %>% left_join(.,ref_file %>% select(-map)) 

a1 = c(sprintf("%s_alt", pops), "talt")
a2 = c(sprintf("%s_ref", pops), "tref")

freqs = lapply(1: (length(pops) + 1), 
               function(i) run_snp[a1[i]] / (run_snp[a1[i]] + run_snp[a2[i]])) %>% 
  bind_cols %>% mutate(random_read = run_snp$random_read) %>%
  as.data.frame

for(pop in pops){
  varname_ = sprintf("%s_alt", pop)
  new_varname_ = sprintf("%s_match", pop)
  random_x = rbinom(length(freqs[,varname_]),1,freqs[,varname_])
  freqs = freqs %>% mutate( !!new_varname_ := ifelse(random_read == random_x, 1, 0))
  
}


run_snp = cbind(run_snp,freqs[,c(sprintf("%s_match", pops))])

for(pop in pops){
  varname_ = sprintf("%s_match", pop)
  new_varname_ = sprintf("n_%s", pop)
  run_snp_n_match_bin = run_snp %>% group_by(bin) %>% dplyr::summarise(!!new_varname_ := sum(!!sym(varname_),na.rm = T))
  bin_called_file = bin_called_file %>% left_join(.,run_snp_n_match_bin[,c("bin",new_varname_)],by=c("id"="bin"))
}


rle_snp = rle_file %>% rowwise %>%
  do(X=get_snp_status(rle_file_x  = ., bin_called_file=bin_called_file,pops=pops)) %>%
  unnest(X) %>% bind_cols(.)

rle_snp$frag_ID <- paste(sample,rle_snp$target,rle_snp$chrom,rle_snp$start,rle_snp$end,sep="_")
rle_snp = rle_snp %>% rename(n_all_snps = n_snp)
bin_called_file <-  bin_called_file %>% left_join(.,rle_snp[,c("frag_ID","n_all_snps",paste0("all_n_",pops))],by="frag_ID")

write.csv(bin_called_file,output_bin,quote = F,row.names = F)

write.csv(rle_snp,output_rle,quote = F,row.names = F)
