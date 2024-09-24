library(tidyverse)

true_path = snakemake@input$true
obs_path = snakemake@input$obs

sample_name = snakemake@wildcards$sample
asc_true = snakemake@wildcards$Tasc
asc_obs = snakemake@wildcards$snpset
downsampling_true = snakemake@wildcards$Tdown
penalty_obs = snakemake@wildcards$penalty
penalty_true = snakemake@wildcards$TruePenalty
downsampling_obs = snakemake@wildcards$downsampled
## only if different maps are specified
map_true = snakemake@wildcards$Tmap
map_obs = snakemake@wildcards$map
type_ = snakemake@params$type
state_ = snakemake@params$state
min_len = snakemake@params$min_len

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
  
  df$target = dplyr::first(overlap$target)
  
  df
}

classify_frags <- function(true,est_path,type_,target_,name,ascertainment,downsampling,penalty,min_len,map_){

  est <- read.csv(est_path)
  est$chrom <- as.character(est$chrom)
  est = est %>% filter(type == type_)  %>% filter(target == target_)  %>% filter(map_len >= min_len)
  
  
  if(nrow(true)> 0){
    # hits records for each true fragment, if there are any overlapping estimated hits
    true_hits = true %>% rowwise %>%
      do(X=get_est_from_true_pos(true_row=., est=est)) %>%
      unnest(X) %>% bind_cols(true) 
    
    
    #' false negatives are all true frags not overlapping anything
    false_negatives = true_hits %>% filter(n_overlap==0) %>%
      mutate(len=pos_len, CLS="FN") %>%
      dplyr::select(len, CLS, start=pos, end=pos_end,map_start=map, map_end=map_end, chrom)
    
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
                target=dplyr::first(target...8), 
                merge_gap=max(pos_end) - min(pos) - sum(pos_end-pos),
                start=min(pos), end=max(pos_end),
                map_start=min(map), map_end=max(map_end),
                Ascertainment = dplyr::first(True_Ascertainment),
                downsampling = dplyr::first(True_downsampling),
                sample = dplyr::first(True_sample)
      ) %>%
      mutate(err_start= est_start - start,
             err_end = est_end - end,
             len=end-start) %>% 
      dplyr::select(chrom, len, start, end, map_start,map_end,est_start, est_end, 
             est_gap, merge_gap, n_overlap, target,
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
  data$Obs_sample <- as.character(name)
  data$Obs_downsampling <-  as.numeric(downsampling)
  data$Obs_Ascertainment <-  as.character(ascertainment)
  data$Obs_penalty <-  as.numeric(penalty)
  data$Obs_map <-  as.character(map_)
  return(data)
}


true <- read.csv(true_path)
true$chrom <- as.character(true$chrom)
true = true %>% filter(type == type_) %>% filter(target == "NEA") %>% filter(map_len >= min_len)
true$True_sample <- paste(sample_name,asc_true,penalty_true,sep="_")
true$True_downsampling <-  as.numeric(downsampling_true)
true$True_Ascertainment <-  asc_true
true$True_penalty <-  as.numeric(penalty_true)
true$True_map <-  as.character(map_true)

obs_res = classify_frags(true = true, 
                                 est_path = obs_path,
                                 type_= type_,
                                 target_ = state_,
                                 name = paste(sample_name,asc_obs,penalty_obs,sep="_"),
                                 ascertainment = asc_obs,
                                 downsampling = as.numeric(downsampling_obs),
                                 penalty = penalty_obs,
                                 min_len = min_len,
                                 map_ = as.character(map_obs))

obs_res$True_sample <- paste(sample_name,asc_true,penalty_true,sep="_")
obs_res$True_downsampling <-  as.numeric(downsampling_true)
obs_res$True_Ascertainment <-  asc_true
obs_res$True_penalty <-  as.numeric(penalty_true)
obs_res$True_map <-  as.character(map_true)

obs_res = obs_res 


write.csv(obs_res,outfile,row.names = F)
