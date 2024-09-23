library(tidyverse)
library(readr)


true_input = snakemake@input$true_rle
ind_ = as.numeric(snakemake@wildcards$id)
rep_ = snakemake@wildcards$rep
rle_admixfrog = snakemake@input$rle_admixfrog
hmmix_rle = snakemake@input$hmmix_hap1matched
min_len = as.numeric(snakemake@wildcards$minlen)
source_ = as.character(snakemake@wildcards$source)
hmmix_matched = as.character(snakemake@wildcards$hmmix_matched)
cutoff_ = as.numeric(snakemake@wildcards$hmmix_cutoff)
outfile = snakemake@output$res_obs


get_est_from_true_pos <- function(true_row, est){
  #       t0------t1
  #    e0----e1
  #            e0----e1
  #   e0---------------e1
  #         e0-e1
  overlap = est %>% filter(pos<true_row$pos_end, 
                           pos_end > true_row$pos,
                           chrom == true_row$chrom
  )
  n_overlap = nrow(overlap)
  suppressWarnings({
    est_start = min(overlap$pos)
    est_end = max(overlap$pos_end)
    est_src = dplyr::first(overlap$target)
    est_map_start =  min(overlap$map)
    est_map_end =  min(overlap$map_end)
    ids = list(overlap$id)
    est_gap = max(overlap$pos_end) - min(overlap$pos) - 
      sum(overlap$pos_end - overlap$pos)
  })
  df = tibble(n_overlap, est_start, est_end, est_gap,est_map_start,est_map_end)
  df[df$n_overlap == 0] = 0
  df$est_ids = ids
  if(is.null(df$est_ids[[1]]))df$est_ids = lapply(1:nrow(df), function(x)tibble(est_ids=integer()))
  
  df$est_target = dplyr::first(overlap$target)
  
  df
}

classify_frags <- function(true,est){
  
  
  
  if(nrow(true)> 0){
    # hits records for each true fragment, if there are any overlapping estimated hits
    true_hits = true %>% rowwise %>%
      do(X=get_est_from_true_pos(true_row=., est=est)) %>%
      unnest(X) %>% bind_cols(true) 
    
    
    #' false negatives are all true frags not overlapping anything
    
    false_negatives = true_hits %>% filter(n_overlap==0)
      
      if(nrow(false_negatives) > 0){
        false_negatives = false_negatives %>%
          mutate(len=pos_len, CLS="FN") %>%
        dplyr::select(len, CLS, start=pos, end=pos_end,map_start=map, map_end=map_end, chrom)
      
    }

    #ids of all called fragments overlapping at least one true frag
    est_hit_ids = true_hits %>% unnest(est_ids) %>% dplyr::select(est_ids) %>% unlist
    
    #' false positives are all the called frags that do not overlap any true frags
    false_positives = est %>% filter(!id %in% est_hit_ids) %>%
      mutate(len=pos_end - pos, CLS="FP") %>%
      dplyr::select(len, CLS, chrom, est_target=target, est_start=pos, est_end=pos_end , est_map_start=map, est_map_end=map_end) 
    
    est_hits = true_hits %>% 
      unnest(est_ids) %>% 
      group_by(est_ids) %>%
      summarize(n_merged=n(),id=list(id), 
                #diploid_id=first(diploid_id), 
                chrom=dplyr::first(chrom), 
                est_start=min(est_start),  est_end=max(est_end),
                est_map_start=min(est_map_start),  est_map_end=max(est_map_end),
                est_gap=sum(est_gap),
                n_overlap=sum(n_overlap),
                #from=dplyr::first(from), 
                est_target=dplyr::first(est_target), 
                source = source,
                merge_gap=max(pos_end) - min(pos) - sum(pos_end-pos),
                start=min(pos), end=max(pos_end),
                map_start=min(map), map_end=max(map_end),
      ) %>%
      mutate(err_start= est_start - start,
             err_end = est_end - end,
             len=end-start) %>% 
      dplyr::select(chrom, len, start, end, map_start,map_end,est_start, est_end, 
                    est_gap, merge_gap, n_overlap, est_target,source,
                    n_merged, err_start, err_end, est_map_start,est_map_end) %>%
      distinct
    
    
    
    if(nrow(est_hits) > 0){
      est_hits$CLS = NA
    } else{
      est_hits$CLS = c()
    }
    
    
    
    H = est_hits
    H$CLS[H$n_overlap * H$n_merged == 1 ] = "TP"
    H$CLS[H$n_overlap > 1] = "GAP"
    H$CLS[H$n_merged > 1] = "MERGED"
    H$CLS[H$n_merged > 1 & H$merge_gap <= 0 ] = "OVERLAP"
    data = bind_rows(H,false_negatives, false_positives) %>%
      mutate(len2=round(len / 1e6, 2) * 1e6)
  } else {
    data = est %>% 
      mutate(len=pos_end - pos, CLS="FP") %>%
      dplyr::select(len, CLS, est_start=pos, est_end=pos_end) 
  }
  
  if('est_ids' %in% colnames(data)){
    data = data %>% ungroup() %>% dplyr::select(!est_ids)
    
  }
  

  return(data)
}

true <- read.csv(true_input) %>%
  filter(sample_name == paste0("target",ind_),source == source_) %>% 
  mutate(chrom = as.character(chrom),map = start*1e-6, map_end = end * 1e-6,pos = start, pos_end = end) %>%
  mutate(map_len = map_end - map,pos_len = pos_end - pos) %>%
  filter(map_len >= min_len) 

admixfrog_rle <- read.csv(rle_admixfrog) %>%
  filter(type == "state",target != "AFR") %>% mutate(chrom = paste0("chr",chrom)) %>%
  mutate(chrom = as.character(chrom),map = pos*1e-6, map_end = pos_end * 1e-6) %>%
  mutate(map_len = map_end - map)

hmmix_rle <- read.csv(hmmix_rle)

if(hmmix_matched == "Matched"){
  hmmix_rle <-  hmmix_rle %>%
    filter(target != "unknown") %>%
    mutate(chrom = as.character(chrom),map = start*1e-6, map_end = end * 1e-6,pos = start, pos_end = end) %>%
    mutate(map_len = map_end - map,pos_len = pos_end - pos) %>%
    rownames_to_column(.,"id") 
} else {
  hmmix_rle <-  hmmix_rle %>%
    mutate(chrom = as.character(chrom),map = start*1e-6, map_end = end * 1e-6,pos = start, pos_end = end) %>%
    mutate(map_len = map_end - map,pos_len = pos_end - pos) %>%
    rownames_to_column(.,"id") 
}





obs_res_admixfrog = classify_frags(true = true, 
                         est = admixfrog_rle)

obs_res_admixfrog <- obs_res_admixfrog %>% 
  mutate(ind = paste0("target",ind_),rep = rep_,method = "admixfrog",min_len = min_len, min_pp_cutoff = NA, matched = NA) 

obs_res_hmmix = classify_frags(true = true, 
                         est = hmmix_rle)

obs_res_hmmix <- obs_res_hmmix %>% 
  mutate(ind = paste0("target",ind_),rep = rep_,method = "hmmix",min_len = min_len, min_pp_cutoff = cutoff_, matched = hmmix_matched) 


obs_res <- rbind(obs_res_admixfrog,obs_res_hmmix) 

write.csv(obs_res,outfile,row.names = F)
