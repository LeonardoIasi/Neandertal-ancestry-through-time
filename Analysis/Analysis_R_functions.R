suppressPackageStartupMessages({
  library("gplots")
  library(reshape2)
  library(dendextend)
  library(paletteer)
  library(tidyverse)
  library(rethinking)
  library(ggrepel)
  library(data.table)
  library(ggpubr)
  library(cowplot)
  library("DEoptim")
  library(VGAM)
  library(DPQ)
  library(parallel)
  library(foreach)
  library(doMC)
  library(binom)
})


'%notin%' <- Negate('%in%')


get_rle_files <- function(directory_path, split_char,type_,target_,chrom_,sep_char="\t",no_filter=F){
  for(i in 1:length(directory_path)){
    if(i == 1){
      smp_name = strsplit(x = tail(unlist(strsplit(directory_path[i],split = "/")),1),split = split_char)[[1]][1]
      print(smp_name)
      if(no_filter){
        print("No filtering")
        All_rle = read.csv(directory_path[i],sep=sep_char) %>%  mutate(sample = smp_name) %>% mutate(chrom = as.character(chrom))
      } else {
        All_rle = read.csv(directory_path[i],sep=sep_char) %>% filter(type %in% type_) %>% 
          filter(target %in% target_) %>% 
          filter(chrom %in% chrom_) %>% mutate(sample = smp_name) %>% mutate(chrom = as.character(chrom))
      }
    } else {
      smp_name = strsplit(x = tail(unlist(strsplit(directory_path[i],split = "/")),1),split = split_char)[[1]][1]
      print(smp_name)
      if(no_filter){
        print("No filtering")
        sample_x = read.csv(directory_path[i],sep=sep_char)  %>% mutate(chrom = as.character(chrom))
      } else {
        sample_x = read.csv(directory_path[i],sep=sep_char)  %>% filter(type %in% type_) %>% 
          filter(target %in% target_) %>% 
          filter(chrom %in% chrom_) %>% mutate(chrom = as.character(chrom))
      }

      if(length(sample_x$id) == 0){
        next
      } else{
        sample_x  = sample_x  %>% mutate(sample = smp_name)
        #All_rle = rbind(All_rle,sample_x, fill = TRUE)
        All_rle = dplyr::bind_rows(All_rle,sample_x)
      }
    }
  }
  return(All_rle)
}

write_bed_from_All_rle <- function(rle_file,type_,target_,chrom_,path,min_len,min_bp_len,min_SNP,optional_col=NULL){
  file_name = c()
  for(smp in unique(rle_file$sample)){
    smp_name = smp
    print(smp_name)
    if(is.null(optional_col)){
      All_rle = rle_file %>% filter(sample == smp) %>% filter(type %in% type_) %>% 
        filter(target %in% target_) %>% 
        filter(chrom %in% chrom_) %>% mutate(sample = smp_name) %>%
        drop_na() %>%
        filter(map_len >= min_len, pos_len >= min_bp_len, n_all_snps >= min_SNP) %>%
        mutate(frag_ID = paste(sample,target,chrom,start,end,sep="_")) %>%
        dplyr::select(chrom,pos,pos_end,sample) %>% mutate(chrom = ifelse(chrom == "X","23",chrom)) %>% 
        mutate_at(c('chrom'), as.numeric) %>% arrange(.,chrom,pos,pos_end)
    } else{
      All_rle = read.csv(directory_path[i]) %>% filter(type %in% type_) %>% 
        filter(target %in% target_) %>% 
        filter(chrom %in% chrom_) %>% mutate(sample = smp_name) %>%
        drop_na() %>%
        filter(map_len >= min_len, pos_len >= min_bp_len, n_all_snps >= min_SNP) %>%
        mutate(frag_ID = paste(sample,target,chrom,start,end,sep="_")) %>%
        dplyr::select(chrom,pos,pos_end,sample,!!optional_col) %>% mutate(chrom = ifelse(chrom == "X","23",chrom)) %>% 
        mutate_at(c('chrom'), as.numeric) %>% arrange(.,chrom,pos,pos_end)
    }
    
    file_name_x = paste(path,"/",smp_name,"_", paste(target_, collapse='_'),"_",paste(type_, collapse='_'),"min_len_",min_len,"_minSNP_",min_SNP,".bed",sep="")  
    write.table(All_rle,file_name_x,quote = F,sep = "\t",row.names = F,col.names = F)
    file_name = c(file_name,file_name_x)
  }
  return(file_name)
}


write_bed_from_rle <- function(directory_path, split_char,type_,target_,chrom_,path,min_len,min_bp_len,min_SNP,optional_col=NULL){
  file_name = c()
  for(i in 1:length(directory_path)){
    smp_name = strsplit(x = tail(unlist(strsplit(directory_path[i],split = "/")),1),split = split_char)[[1]][1]
    print(smp_name)
    if(is.null(optional_col)){
      All_rle = read.csv(directory_path[i]) %>% filter(type %in% type_) %>% 
        filter(target %in% target_) %>% 
        filter(chrom %in% chrom_) %>% mutate(sample = smp_name) %>%
        drop_na() %>%
        filter(map_len >= min_len, pos_len >= min_bp_len, n_all_snps >= min_SNP) %>%
        mutate(frag_ID = paste(sample,target,chrom,start,end,sep="_")) %>%
        dplyr::select(chrom,pos,pos_end,sample) %>% mutate(chrom = ifelse(chrom == "X","23",chrom)) %>% 
        mutate_at(c('chrom'), as.numeric) %>% arrange(.,chrom,pos,pos_end)
    } else{
      All_rle = read.csv(directory_path[i]) %>% filter(type %in% type_) %>% 
        filter(target %in% target_) %>% 
        filter(chrom %in% chrom_) %>% mutate(sample = smp_name) %>%
        drop_na() %>%
        filter(map_len >= min_len, pos_len >= min_bp_len, n_all_snps >= min_SNP) %>%
        mutate(frag_ID = paste(sample,target,chrom,start,end,sep="_")) %>%
        dplyr::select(chrom,pos,pos_end,sample,!!optional_col) %>% mutate(chrom = ifelse(chrom == "X","23",chrom)) %>% 
        mutate_at(c('chrom'), as.numeric) %>% arrange(.,chrom,pos,pos_end)
    }

    file_name_x = paste(path,"/",smp_name,"_", paste(target_, collapse='_'),"_",paste(type_, collapse='_'),"min_len_",min_len,"_minSNP_",min_SNP,".bed",sep="")  
    write.table(All_rle,file_name_x,quote = F,sep = "\t",row.names = F,col.names = F)
    file_name = c(file_name,file_name_x)
  }
  return(file_name)
}

merge_bins_fn <- function(dir_path_anc,dir_path_pd,
                          called_Ancestry,map,
                          bin_called_min_len_anc,
                          bin_called_min_len_pd,
                          bin_called_min_len_pos_anc,
                          bin_called_min_len_pos_pd,min_n_all_SNP_anc,min_n_all_SNP_pd,
                          min_n_SNP_anc,min_n_SNP_pd,out_dir = NULL) {
  anc_x <- read.csv(paste0(dir_path_anc,"/",map,"/","All_bins_",called_Ancestry,"_CurratedEMHPublished_archaicadmixtureAPX.bin_called0.25_min_len",bin_called_min_len_anc,"_min_len_pos",bin_called_min_len_pos_anc,"_min_n_all_SNP",min_n_all_SNP_anc,"_min_n_SNP",min_n_SNP_anc,".xz"))
  pd_x <- read.csv(paste0(dir_path_pd,"/",map,"/","All_bins_",called_Ancestry,"_PresentDaySGDP_archaicadmixtureAPX.bin_called0.25_min_len",bin_called_min_len_pd,"_min_len_pos",bin_called_min_len_pos_pd,"_min_n_all_SNP",min_n_all_SNP_pd,"_min_n_SNP",min_n_SNP_pd,".xz"))
  res_x <- anc_x %>% dplyr::select(!c(chrom,pos)) %>% inner_join(pd_x,.,by=c("id","map"))
  if(!is.null(out_dir)){
    write.csv(res_x,
              paste0(out_dir,"All_merged_bins_",called_Ancestry,"archaicadmixtureAPX.bin_called_map_",map,"_penatly_0.25_min_len",bin_called_min_len_anc,"_",bin_called_min_len_pd,"_min_len_pos",bin_called_min_len_pos_anc,"_",bin_called_min_len_pos_pd,"_min_n_all_SNP",min_n_all_SNP_anc,"_",min_n_all_SNP_pd,"_min_n_SNP",min_n_SNP_anc,"_",min_n_SNP_pd,".xz"),
              quote = F,row.names = F)
    
  } else {
    return(res_x)
  }
  
}

merge_rle_fn <- function(dir_path_anc,dir_path_pd,
                         called_Ancestry,map,
                         called_type,
                         penalty,
                         chrom_,
                         bin_called_min_len_anc,
                         bin_called_min_len_pd,
                         bin_called_min_len_pos_anc,
                         bin_called_min_len_pos_pd,
                         min_n_all_SNP_anc,min_n_all_SNP_pd,
                         split_char,
                         out_dir = NULL,
                         no_filter = F) {
  dir_path_EMH <- list.files(path=paste0(dir_path_anc,map),pattern = paste0("rle",penalty,".xz.anno_all_types$"),full.names = T)
  
  
  dir_path_SGDP <- list.files(path=paste0(dir_path_pd,map),pattern = paste0("rle",penalty,".xz.anno_all_types$"),full.names = T)
  if(no_filter){
    EMH_rle_all <- get_rle_files(dir_path_EMH,split_char,type_ = called_type,target_ = called_Ancestry,chrom_,sep_char = ",",no_filter = no_filter) 
    SGDP_rle_all <- get_rle_files(dir_path_SGDP,split_char,type_ = called_type,target_ = called_Ancestry,chrom_,sep_char=",",no_filter = no_filter) %>% 
      mutate(sample = gsub("-", "\\.", sample))
  } else {
    EMH_rle_all <- get_rle_files(dir_path_EMH,split_char,type_ = called_type,target_ = called_Ancestry,chrom_,sep_char = ",",no_filter = no_filter) %>% 
      drop_na() %>% filter(map_len >= bin_called_min_len_anc, pos_len >= bin_called_min_len_pos_anc,n_all_snps >= min_n_all_SNP_anc)
    SGDP_rle_all <- get_rle_files(dir_path_SGDP,split_char,type_ = called_type,target_ = called_Ancestry,chrom_,sep_char=",",no_filter = no_filter) %>% 
      filter(map_len >= bin_called_min_len_pd, pos_len >= bin_called_min_len_pos_pd,n_all_snps >= min_n_all_SNP_pd) %>%
      mutate(sample = gsub("-", "\\.", sample))
  }
  
  
  All_rle <- rbind(SGDP_rle_all,EMH_rle_all)
  if(!is.null(out_dir)){
    called_Ancestry <- paste(called_Ancestry,collapse = "_")
    write.csv(All_rle,
              paste0(out_dir,"All_merged_rle_",called_Ancestry,"archaicadmixtureAPX.bin_called_map_",map,"_penatly_0.25_min_len",bin_called_min_len_anc,"_",bin_called_min_len_pd,"_min_len_pos",bin_called_min_len_pos_anc,"_",bin_called_min_len_pos_pd,"_min_n_all_SNP",min_n_all_SNP_anc,"_",min_n_all_SNP_pd,".xz"),
              quote = F,row.names = F)
    
  } else {
    return(All_rle)
  }
  
}

group_overlaps_fn <- function(fragments_data){
  DT_xx <- as.data.table(fragments_data)
  
  i = 1
  for(chr in unique(fragments_data$chrom)){
    
    print(paste("Processing chromosome: ",chr))
    DT_xx_chr = DT_xx %>% filter(chrom == chr)
    fragments_data_chr = fragments_data %>% filter(chrom == chr)
    
    DT_xx_chr.int <- as.data.table(
      intervals::interval_union( 
        intervals::Intervals( as.matrix( DT_xx_chr[, c("pos","pos_end")] ) ) , 
        check_valid = TRUE ) )
    
    colnames(DT_xx_chr.int) <- c("start", "end" )
    DT_xx_chr.int$group_id <- paste(chr,rownames(DT_xx_chr.int),sep = "_")
    
    setkey( DT_xx_chr, pos, pos_end)
    setkey( DT_xx_chr.int, start, end)
    fragments_data_chr <- foverlaps( DT_xx_chr.int, DT_xx_chr )
    group_counts_chr = fragments_data_chr %>% dplyr::count(group_id)
    
    fragments_data_chr = inner_join(fragments_data_chr,group_counts_chr,by="group_id") %>%
      mutate(group_id_len = end - start)
    fragments_data_chr$overlapp_density <- fragments_data_chr$n/(fragments_data_chr$end - fragments_data_chr$start)
    
    if(i == 1){
      fragments_data_res = fragments_data_chr
      
    }else{
      fragments_data_res = rbind(fragments_data_res,fragments_data_chr)
    }
    
    i = i + 1
    
    
    
  }
  
  
  
  return(fragments_data_res)
  
}

merge_consecutive_bins_fn <- function(data){
  DT_xx <- as.data.table(data)
  
  i = 1
  for(chr in unique(data$chrom)){
    
    print(paste("Processing chromosome: ",chr))
    DT_xx_chr = DT_xx %>% filter(chrom == chr)
    data_chr = data %>% filter(chrom == chr)
    
    DT_xx_chr.int <- as.data.table(
      intervals::interval_union( 
        intervals::Intervals( as.matrix( DT_xx_chr[, c("start","end")] ) ) , 
        check_valid = TRUE ) )
    
    colnames(DT_xx_chr.int) <- c("region_start", "region_end" )
    DT_xx_chr.int$group_id <- paste(chr,rownames(DT_xx_chr.int),sep = "_")
    
    setkey( DT_xx_chr, start, end)
    setkey( DT_xx_chr.int, region_start, region_end)
    data_chr <- foverlaps( DT_xx_chr.int, DT_xx_chr )
    group_counts_chr = data_chr %>% dplyr::count(group_id)
    
    data_chr = inner_join(data_chr,group_counts_chr,by="group_id") %>%
      mutate(region_len = region_end - region_start) %>% 
      #mutate(region_len = region_end - region_start) %>% dplyr::select(!c(start,end,group_id)) %>%
      distinct()
    
    if(i == 1){
      data_res = data_chr
      
    }else{
      data_res = rbind(data_res,data_chr)
    }
    
    i = i + 1
    
    
    
  }
  
  
  
  return(data_res)
  
}

generat_Superpop_col <- function(factor_sequence,col_pallet){
  population_colors <- unique(cluster_color)
  population_colors <- population_colors[match(factor_sequence,population_colors$clusterF3D),]
  population_colors$clusterF3D <- as.numeric(as.factor(population_colors$clusterF3D))
  return(population_colors$cluster_color)
}

merge_pop_vectors_fn <- function(pop.list,data.matrix,Meta.info,BooleanCoding=T){
  joint_by_pop = data.frame(matrix(NA,nrow = length(data.matrix$chrom),ncol = length(pop.list)))
  i = 0
  for(pop_ in pop.list){
    i = i + 1
    print(paste("merging ",pop_))
    pop_samples = Meta.info %>% filter(pop == pop_) %>% dplyr::select(sample_name)
    joint_by_pop[,i] <- rowMeans(data.matrix[pop_samples$sample_name])
    if(BooleanCoding == T){
      joint_by_pop[,i] <- ifelse(joint_by_pop[,i] > 0,1,0)
    }
    colnames(joint_by_pop)[i] <- pop_ 
  }
  return(joint_by_pop)
}

get_overlaps <- function(x,df,self_overlap=T){
  overlaping_segs_x = df %>% filter(chrom == x$chrom, pos < x$pos_end, pos_end > x$pos)
  sample_names_x = paste(sort(overlaping_segs_x$sample),collapse = ",")
  overlap_id_x <- paste(sort(overlaping_segs_x$frag_ID),collapse = ",")
  count_x <- length(overlaping_segs_x$sample)
  if(self_overlap == F){
    
    sample_names_x = paste(sort(overlaping_segs_x$sample[overlaping_segs_x$sample != x$sample]),collapse = ",")
    overlap_id_x <- paste(sort(overlaping_segs_x$frag_ID[overlaping_segs_x$frag_ID != x$frag_ID]),collapse = ",")
    count_x <- length(overlaping_segs_x$sample[overlaping_segs_x$sample != x$sample])
  }
  x$overlap <-  sample_names_x
  x$overlap_frag_ID <-  overlap_id_x
  x$overlap_count <- count_x
  return(x)
  
}

get_all_segment_overlaps <- function(data,self_overlap=T){
  i = 1
  col_names <- c(colnames(data),"overlap_samples","overlap_frag_ID","overlap_count")
  for(chr in unique(data$chrom)){
    print(paste("Processing chromosome: ",chr))
    data_red = data %>% filter(chrom == chr)
    
    res_overlap_list = data_red %>% rowwise %>%
      do(X=get_overlaps(x=., df=data_red,self_overlap))
    if(i == 1){
      res_overlap =  data.frame(matrix(unlist(res_overlap_list), ncol = max(lengths(res_overlap_list$X)), byrow = T))
      colnames(res_overlap) <- col_names
    }else{
      res_overlap_x =  data.frame(matrix(unlist(res_overlap_list), ncol = max(lengths(res_overlap_list$X)), byrow = T))
      colnames(res_overlap_x) <- col_names
      res_overlap = rbind(res_overlap,res_overlap_x)
    }
    i = i + 1 
  }
  return(res_overlap)
}

## Nea seg cor vs 1240 fn

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

reorder_F3 <- function(F3){
  # calculate distance from F3s
  dd <- dist(F3)
  hc <- hclust(dd)
  F3 <-F3[hc$order, hc$order]
}

annotate_matrix <- function(long_data_matrix,value_name,cluster_used,match_to_matrix = NULL){
  
  long_data_matrix_anno <- long_data_matrix %>% inner_join(.,Joint_Meta[,c("ID",cluster_used,"SuperclusterF3D","ML_BP_Mean")],by=c("Var1"="ID")) %>% dplyr::rename(Var1_cluster = !!cluster_used,Var1_ML_BP_Mean = ML_BP_Mean) %>%
    inner_join(.,Joint_Meta[,c("ID",cluster_used,"SuperclusterF3D","ML_BP_Mean")],by=c("Var2"="ID")) %>% dplyr::rename(Var2_cluster = !!cluster_used,Var2_ML_BP_Mean = ML_BP_Mean) %>% distinct() %>%
    filter(Var1_cluster != "Africa" & Var2_cluster != "Africa") 
  
  if(is.null(match_to_matrix)) {
    match_to_matrix = long_data_matrix
  }
  
  long_data_matrix_anno = long_data_matrix_anno[ order(match(long_data_matrix_anno$Var1, match_to_matrix$Var1)), ]
  long_data_matrix_anno$Var1 <- ordered(long_data_matrix_anno$Var1, levels = levels(match_to_matrix$Var1))
  long_data_matrix_anno = long_data_matrix_anno[ order(match(long_data_matrix_anno$Var2, match_to_matrix$Var2)), ]
  long_data_matrix_anno$Var2 <- ordered(long_data_matrix_anno$Var2, levels = levels(match_to_matrix$Var2))
  colnames(long_data_matrix_anno)[colnames(long_data_matrix_anno) == "value"] <- value_name
  
  return(long_data_matrix_anno)
}

normalize_standardize_values <- function(long_data_matrix_anno,value_name){ 
  vale_wo_offsprings <- long_data_matrix_anno %>% filter(Var1 != "AfanasievoSon1" & Var2 != "AfanasievoSon1" & Var1 != "AfanasievoSon2" & Var2 != "AfanasievoSon2") 
  
  long_data_matrix_anno[,paste(value_name,"normalized",sep="_")] <- (long_data_matrix_anno[,value_name] - min(vale_wo_offsprings[,value_name],na.rm = T)) / (max(vale_wo_offsprings[,value_name],na.rm = T) - min(long_data_matrix_anno[,value_name],na.rm = T) )
  
  long_data_matrix_anno[,paste(value_name,"normalized",sep="_")]  <- ifelse(long_data_matrix_anno[,paste(value_name,"normalized",sep="_")] < 0, 0, ifelse(long_data_matrix_anno[,paste(value_name,"normalized",sep="_")] > 1,1 ,long_data_matrix_anno[,paste(value_name,"normalized",sep="_")]))
  
  long_data_matrix_anno[,paste(value_name,"standardize",sep="_")] <- (long_data_matrix_anno[,value_name] - mean(vale_wo_offsprings[,value_name],na.rm = T)) / sd(vale_wo_offsprings[,value_name],na.rm = T)
  
  return(long_data_matrix_anno)
}

ggHeatmap <- function(annotated_matrix,min_limit,max_limit,col_value_name){
  col_value_name <- sym(col_value_name)
  ggplot(data = annotated_matrix, aes(x=Var1, y=Var2, fill=!!col_value_name)) + 
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(min_limit,max_limit), space = "Lab")+
    geom_point(data= annotated_matrix, aes(x=Var1, y=Var2, col=Var1_cluster),alpha=0) +
    xlab("") + 
    ylab("") +
    scale_color_manual(name='Genetic Clusters',
                       labels = cluster_color$clusterF3D,values = cluster_color$cluster_color) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))
           #,fill = 'none'
    )+
    theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                     size = 8, hjust = 1, colour = generat_Superpop_col(annotated_matrix$Var1_cluster[!duplicated(annotated_matrix$Var1)],cluster_color)),
          axis.text.y = element_text( vjust = 1, 
                                      size = 8, hjust = 1, colour = generat_Superpop_col(annotated_matrix$Var2_cluster[!duplicated(annotated_matrix$Var2)],cluster_color)))
  
}

## Unique segment
Reshape_pairwise.t.test_PopClust_fn <- function(pairwise.t.test_res,n_row=12){
  # this function gets a matrix from pairwise.t.test and assumes the population clusters are alphabetically ordered (does not work for a differnt set of pop clusters)
  pairwise.t.test_res_pvals_row1 <- matrix(data = c(1,rep(NA,n_row)),nrow = 1)
  colnames(pairwise.t.test_res_pvals_row1) <- c(colnames(pairwise.t.test_res$p.value))
  rownames(pairwise.t.test_res_pvals_row1) <- "AncientEAS"
  pairwise.t.test_res_pvals_colk <- matrix(data = c(rep(NA,n_row+1),1),ncol = 1)
  colnames(pairwise.t.test_res_pvals_colk) <- c("Yamnaya")
  rownames(pairwise.t.test_res_pvals_colk) <- c("AncientEAS",rownames(pairwise.t.test_res$p.value))
  pairwise.t.test_res_pvals <- cbind(rbind(pairwise.t.test_res_pvals_row1,pairwise.t.test_res$p.value),pairwise.t.test_res_pvals_colk)
  diag(pairwise.t.test_res_pvals) <- 1
  for(i in 1:nrow(pairwise.t.test_res_pvals)){
    for(j in 1:ncol(pairwise.t.test_res_pvals)){
      if(is.na(pairwise.t.test_res_pvals[i,j])) pairwise.t.test_res_pvals[i,j] <- pairwise.t.test_res_pvals[j,i]
    }
  }
  return(pairwise.t.test_res_pvals)
}

## Matching fn

match_fn <- function(MatchingRates_file_anno){
  MatchingRates_file_anno_class <- MatchingRates_file_anno %>% filter(target == "NEA") %>% filter(chrom != "X", n_all_snps > 10) %>%
    filter(!if_all(all_of(c("prop_pHcount_matching_Vindija33.19","prop_pHcount_matching_Altai","prop_pHcount_matching_Chagyrskaya","prop_pHcount_matching_Denisova")), is.na))
  
  Highest_matching_archaic = MatchingRates_file_anno_class %>% dplyr::select(prop_pHcount_matching_Vindija33.19,prop_pHcount_matching_Altai,prop_pHcount_matching_Chagyrskaya,prop_pHcount_matching_Denisova) %>%
    replace(is.na(.), 0) %>% dplyr::rename(Vindija33.19 = prop_pHcount_matching_Vindija33.19,
                                           Altai = prop_pHcount_matching_Altai,
                                           Chagyrskaya = prop_pHcount_matching_Chagyrskaya,
                                           Denisova = prop_pHcount_matching_Denisova) %>% apply(.,1,function(x) paste(sort(names(which(x==max(x)))),collapse = ", "))
  
  Highest_matching_archaic_value = MatchingRates_file_anno_class %>% dplyr::select(prop_pHcount_matching_Vindija33.19,prop_pHcount_matching_Altai,prop_pHcount_matching_Chagyrskaya,prop_pHcount_matching_Denisova) %>%
    replace(is.na(.), 0) %>%  mutate(Highest_matching_archaic_value = pmax(prop_pHcount_matching_Vindija33.19,prop_pHcount_matching_Altai,prop_pHcount_matching_Chagyrskaya,prop_pHcount_matching_Denisova))
  
  MatchingRates_file_anno_class$Highest_matching_archaic <- Highest_matching_archaic
  MatchingRates_file_anno_class$Highest_matching_archaic_value <- Highest_matching_archaic_value$Highest_matching_archaic_value
  
  MatchingRates_file_anno_class <- MatchingRates_file_anno_class %>% mutate(Highest_matching_archaic = ifelse(Highest_matching_archaic_value == 0,NA,Highest_matching_archaic))
  
  MatchingRates_file_anno_class = MatchingRates_file_anno_class %>% mutate(tree = ifelse((Highest_matching_archaic == "Altai, Chagyrskaya, Vindija33.19" | Highest_matching_archaic == "Chagyrskaya, Vindija33.19" | Highest_matching_archaic == "Altai, Chagyrskaya, Denisova, Vindija33.19" | 
                                                                                            Highest_matching_archaic == "Altai" |
                                                                                            Highest_matching_archaic == "Chagyrskaya" |
                                                                                            Highest_matching_archaic == "Denisova" |
                                                                                            Highest_matching_archaic == "Vindija33.19" ),"cordant","discordant"))
  return(MatchingRates_file_anno_class)
}

match_BT_fn <- function(MatchingRates_file_anno){
  MatchingRates_file_anno_class <- MatchingRates_file_anno %>% filter(target == "NEA") %>% filter(chrom != "X", n_all_snps > 10) %>%
    filter(!if_all(all_of(c("BT_prop_pHcount_matching_Vindija33.19_mean","BT_prop_pHcount_matching_Altai_mean","BT_prop_pHcount_matching_Chagyrskaya_mean","BT_prop_pHcount_matching_Denisova_mean")), is.na))
  
  Highest_matching_archaic = MatchingRates_file_anno_class %>% dplyr::select(BT_prop_pHcount_matching_Vindija33.19_mean,BT_prop_pHcount_matching_Altai_mean,BT_prop_pHcount_matching_Chagyrskaya_mean,BT_prop_pHcount_matching_Denisova_mean) %>%
    replace(is.na(.), 0) %>% dplyr::rename(Vindija33.19 = BT_prop_pHcount_matching_Vindija33.19_mean,
                                           Altai = BT_prop_pHcount_matching_Altai_mean,
                                           Chagyrskaya = BT_prop_pHcount_matching_Chagyrskaya_mean,
                                           Denisova = BT_prop_pHcount_matching_Denisova_mean) %>% apply(.,1,function(x) paste(sort(names(which(x==max(x)))),collapse = ", "))
  
  Highest_matching_archaic_value = MatchingRates_file_anno_class %>% dplyr::select(BT_prop_pHcount_matching_Vindija33.19_mean,BT_prop_pHcount_matching_Altai_mean,BT_prop_pHcount_matching_Chagyrskaya_mean,BT_prop_pHcount_matching_Denisova_mean) %>%
    replace(is.na(.), 0) %>%  mutate(Highest_matching_archaic_value = pmax(BT_prop_pHcount_matching_Vindija33.19_mean,BT_prop_pHcount_matching_Altai_mean,BT_prop_pHcount_matching_Chagyrskaya_mean,BT_prop_pHcount_matching_Denisova_mean))
  
  MatchingRates_file_anno_class$Highest_matching_archaic <- Highest_matching_archaic
  MatchingRates_file_anno_class$Highest_matching_archaic_value <- Highest_matching_archaic_value$Highest_matching_archaic_value
  
  MatchingRates_file_anno_class <- MatchingRates_file_anno_class %>% mutate(Highest_matching_archaic = ifelse(Highest_matching_archaic_value == 0,NA,Highest_matching_archaic))
  
  MatchingRates_file_anno_class = MatchingRates_file_anno_class %>% mutate(tree = ifelse((Highest_matching_archaic == "Altai, Chagyrskaya, Vindija33.19" | Highest_matching_archaic == "Chagyrskaya, Vindija33.19" | Highest_matching_archaic == "Altai, Chagyrskaya, Denisova, Vindija33.19" | 
                                                                                            Highest_matching_archaic == "Altai" |
                                                                                            Highest_matching_archaic == "Chagyrskaya" |
                                                                                            Highest_matching_archaic == "Denisova" |
                                                                                            Highest_matching_archaic == "Vindija33.19" ),"cordant","discordant"))
  return(MatchingRates_file_anno_class)
}

classify_seg_divergence <- function(data){
  data_class <- data  %>%
    filter(!if_all(all_of(c("prop_matching_Vindija33.19","prop_matching_Altai","prop_matching_Chagyrskaya","prop_matching_Denisova")), is.na))
  
  Highest_matching_archaic_reads = data_class %>% dplyr::select(prop_matching_Vindija33.19,prop_matching_Altai,prop_matching_Chagyrskaya,prop_matching_Denisova) %>%
    replace(is.na(.), 0) %>% dplyr::rename(Vindija33.19 = prop_matching_Vindija33.19,
                                           Altai = prop_matching_Altai,
                                           Chagyrskaya = prop_matching_Chagyrskaya,
                                           Denisova = prop_matching_Denisova) %>% apply(.,1,function(x) paste(sort(names(which(x==max(x)))),collapse = ", "))
  
  Highest_matching_archaic_value_reads = data_class %>% dplyr::select(prop_matching_Vindija33.19,prop_matching_Altai,prop_matching_Chagyrskaya,prop_matching_Denisova) %>%
    replace(is.na(.), 0) %>%  mutate(Highest_matching_archaic_value_reads = pmax(prop_matching_Vindija33.19,prop_matching_Altai,prop_matching_Chagyrskaya,prop_matching_Denisova))
  
  Highest_matching_archaic_pH = data_class %>% dplyr::select(prop_pHcount_matching_Vindija33.19,prop_pHcount_matching_Altai,prop_pHcount_matching_Chagyrskaya,prop_pHcount_matching_Denisova) %>%
    replace(is.na(.), 0) %>% dplyr::rename(Vindija33.19 = prop_pHcount_matching_Vindija33.19,
                                           Altai = prop_pHcount_matching_Altai,
                                           Chagyrskaya = prop_pHcount_matching_Chagyrskaya,
                                           Denisova = prop_pHcount_matching_Denisova) %>% apply(.,1,function(x) paste(sort(names(which(x==max(x)))),collapse = ", "))
  
  Highest_matching_archaic_value_pH = data_class %>% dplyr::select(prop_pHcount_matching_Vindija33.19,prop_pHcount_matching_Altai,prop_pHcount_matching_Chagyrskaya,prop_pHcount_matching_Denisova) %>%
    replace(is.na(.), 0) %>%  mutate(Highest_matching_archaic_value_pH = pmax(prop_pHcount_matching_Vindija33.19,prop_pHcount_matching_Altai,prop_pHcount_matching_Chagyrskaya,prop_pHcount_matching_Denisova))
  
  data_class$Highest_matching_archaic_reads <- Highest_matching_archaic_reads
  data_class$Highest_matching_archaic_value_reads <- Highest_matching_archaic_value_reads$Highest_matching_archaic_value_reads
  data_class$Highest_matching_archaic_pH <- Highest_matching_archaic_pH
  data_class$Highest_matching_archaic_value_pH <- Highest_matching_archaic_value_pH$Highest_matching_archaic_value_pH
  
  data_class <- data_class %>% mutate(Highest_matching_archaic_reads = ifelse(Highest_matching_archaic_value_reads == 0,NA,Highest_matching_archaic_reads),
                                      Highest_matching_archaic_pH = ifelse(Highest_matching_archaic_value_pH == 0,NA,Highest_matching_archaic_pH))
  
  return(data_class)
  
}



## Dating fn


## ACov fitting fn

fnx_SP <- function(A,tm,C,offset,dist){
  res_curve <- c()
  for(i in 1:length(offset)){
    res_y = (A[i]*exp(-(dist[[i]])*(tm-offset[i]))+C[i])
    res_curve = c(res_curve,res_y)
  }
  return(res_curve)
}

fnx_SP_singlefit <- function(A,tm,C,offset,dist){
  res_curve = (A*exp(-(dist)*(tm-offset))+C)
  return(res_curve)
}


fnx_TP_singlefit <- function(A1,A2,tm1,tm2,C,offset,dist){
  res_curve = C + (A1*exp(-(dist)*(tm1-offset)) + A2*exp(-(dist)*(tm2-offset)))
  return(res_curve)
}

Joint_Simple_Pulse_ACov_weighted_fn <- function(data.in,minFac,maxIter,JointInterceptsFit){
  
  
  A_est=list(rep(0,length(data.in$offset)))
  tm_est=list(0)
  C_est=list(rep(0,length(data.in$offset)))
  
  A_est[2]=list(data.in$A_init)
  tm_est[2]=list(data.in$tm_init)
  C_est[2]=list(data.in$C_init)
  
  RSS_all <- list(1000,999)
  
  k <- 2
  
  if(JointInterceptsFit == T){
    print("Joint fitting of the amplitude & X-intercept")
  } else {
    print("Independent fit of the amplitude & X-intercept")
  }
  
  #while(abs(RSS_all[[k]] - RSS_all[[k-1]]) > minFac & k <= maxIter + 1) {
  while((RSS_all[[k-1]] - RSS_all[[k]]) > minFac & k <= maxIter + 1) {
    print(paste("iteration: ",k-1))
    
    if(JointInterceptsFit == T){
      
      fit1_Amplitude <- nls(y ~ fnx_SP(A,tm=tm_est[[k]],C=C_est[[k]],offset,dist), start=list(A=A_est[[k]]), algorithm="port",
                            lower=c(rep(1e-6,length(offset))), upper=c(rep(1,length(offset))),
                            control=list(warnOnly=TRUE),data = data.in)
      A_est[k+1] = list(coef(fit1_Amplitude))
      
    } else {
      
      fit1_Amplitude <- c()
      for(a in 1:length(A_est[[k]])){
        y = data.in$ACov_res_list[[a]]$wcov
        x = data.in$ACov_res_list[[a]]$dist
        ts = data.in$offset[a]
        fit1_Amplitude[a] <- coef(nls(y ~ fnx_SP_singlefit(A,tm=tm_est[[k]],C=C_est[[k]][a],offset=ts,dist=x),start=list(A=A_est[[k]][a]), algorithm="port",
                                      lower=c(1e-6), upper=c(1),
                                      control=list(warnOnly=TRUE)))
      }
      A_est[k+1] = list(fit1_Amplitude)
    }
    
    
    fit1_tm <- nls(y ~ fnx_SP(A=A_est[[k+1]],tm,C=C_est[[k]],offset,dist), start=list(tm=tm_est[[k]]), algorithm="port",
                   lower=c(1), upper=c(5000),
                   control=list(warnOnly=TRUE),data = data.in)
    tm_est[k+1] = list(coef(fit1_tm))
    
    
    
    if(JointInterceptsFit == T){
      
      fit1_C <- nls(y ~ fnx_SP(A=A_est[[k+1]],tm=tm_est[[k+1]],C,offset,dist), start=list(C=C_est[[k]]), algorithm="port",
                    lower=c(rep(0,length(offset))), upper=c(rep(1,length(offset))),
                    control=list(warnOnly=TRUE),data = data.in)
      C_est[k+1] = list(coef(fit1_C))
      
    } else {
      
      fit1_C <- c()
      for(c in 1:length(C_est[[k]])){
        y = data.in$ACov_res_list[[c]]$wcov
        x = data.in$ACov_res_list[[c]]$dist
        ts = data.in$offset[c]
        fit1_C[c] <- coef(nls(y ~ fnx_SP_singlefit(A=A_est[[k+1]][c],tm=tm_est[[k+1]],C,offset=ts,dist=x),start=list(C=C_est[[k]][c]), algorithm="port",
                              lower=c(-1), upper=c(1),
                              control=list(warnOnly=TRUE)))
      }
      C_est[k+1] = list(fit1_C)
    }
    
    RSS_all[k+1] = sum(((data.in$y - fnx_SP(A_est[[k+1]],tm_est[[k+1]],C_est[[k+1]],data.in$offset,data.in$dist))^2)* data.in$weight)
    
    k = k + 1
    #print(paste(A_est[[k]],tm_est[[k]],C_est[[k]],RSS_all[[k]]))
  }
  
  return(list(A=A_est[[k]],tm=tm_est[[k]],C=C_est[[k]],RSS=RSS_all[[k]]))
}

Simple_Single_Pulse_ALD_fn <- function(Data,tm_lower,tm_upper){
  xx=Data
  
  dist=xx[,1]
  wcorr=xx[,2]
  
  fm1_exp <- function(x) x[1]*exp(-(dist/100)*x[2])+x[3]
  fm2_exp <- function(x) sum((wcorr-fm1_exp(x))^2)
  fm3_exp <- DEoptim(fm2_exp, lower=c(1e-6,tm_lower,-1), upper=c(1, tm_upper,1), control=list(trace=FALSE))
  
  par1_exp <- fm3_exp$optim$bestmem
  names(par1_exp) <- c("A", "tm","c")
  
  fit1_exp <- nls(wcorr ~ (A*exp(-(dist/100)*tm)+c), start=par1_exp, control=list(maxiter=10000, warnOnly=TRUE)) 
  f_predict = data.frame(dist_M = xx[,1]/100, ALD = predict(fit1_exp))
  
  return(list(fit1_exp,data.frame(tm=coef(fit1_exp)[2], RSS = sum(resid(fit1_exp)^2),lower_trunc=min(xx[,1]),upper_trunc=max(xx[,1])),xx,f_predict))
}


fnx_EP <- function(A,td,tm,C,offset,dist,td_invers=T){
  res_curve <- c()
  for(i in 1:length(offset)){
    if(td_invers){
      td = 1/td  
    }
    tm_corrected = tm - offset[i]
    k = 1/((td/(4*tm_corrected))^2)
    #res_y = A[i] * ( (1 + (tm_corrected/k) * dist[[i]])^(-k) ) + C[i]
    k = 1/k
    res_y = C[i]+A[i]*(1/(1  + ((tm_corrected/ (1/k)) *dist[[i]]) ))^((1/k))
    res_curve = c(res_curve,res_y)
  }
  return(res_curve)
}


fnx_EP_singlefit <- function(A,td,tm,C,offset,dist,td_invers=T){
  if(td_invers){
    td = 1/td  
  }
  tm_corrected = tm - offset
  k = 1/((td/(4*(tm_corrected)))^2)
  #res_curve = A * ( (1 + ((tm_corrected)/k) * dist)^(-k) ) + C
  k = 1/k
  res_curve = C+A*(1/(1  + ((tm_corrected/ (1/k)) *(dist)) ))^((1/k))
  return(res_curve)
}

## fixed tm and td values

GriD_EP_ACov_amplitude_x_inter_fit_weighted <- function(data.in,td_fixed,tm_fixed,A_init,C_init){
  
  fit1_intercepts <- list()
  A_est = c()
  C_est = c()
  for(a in 1:length(A_init)){
    y = data.in$ACov_res_list[[a]]$wcov
    x = data.in$ACov_res_list[[a]]$dist
    ts = data.in$offset[a]
    fit1_intercepts[[a]] <- nls(y ~ fnx_EP_singlefit(A,td=td_fixed,tm=tm_fixed,C,offset=ts,dist=x,td_invers = F),start=list(A=A_init[a],C=C_init[a]), algorithm="port",
                                lower=c(1e-6,0), upper=c(1,1),
                                control=list(warnOnly=TRUE))
    A_est <- c(A_est,coef(fit1_intercepts[[a]])[1])
    C_est <- c(C_est,coef(fit1_intercepts[[a]])[2])
    
    
  }
  
  
  RSS_all = sum(((data.in$y - fnx_EP(A_est,td_fixed,tm_fixed,C_est,data.in$offset,data.in$dist,td_invers = F))^2)  * data.in$weight)
  return(list(A=A_est,td=td_fixed,tm=tm_fixed,C=C_est,RSS=RSS_all))
}


## Segment fitting fn


Two_pulses_EM_SS <- function(data,alpha1_hat,alpha2_hat,t1_hat,t2_hat,num_iterations=100){
  
  for (iteration in 1:num_iterations) {
    # E-step: Calculate the probabilities given the parameters
    p1 <- (alpha1_hat * dexp(data, t1_hat)) / ((alpha1_hat * dexp(data, t1_hat )) + (alpha2_hat * dexp(data, t2_hat)))
    p2 <- 1 - p1
    
    # M-step: Update parameters
    alpha1_hat <- mean(p1)
    alpha2_hat <- mean(p2)
    t1_hat <- 1/ (sum(p1 * data) / sum(p1))
    t2_hat <- 1/ (sum(p2 * data) / sum(p2))
    
    # Print the parameter estimates for this iteration
    #cat("Iteration", iteration, ": alpha1 =", alpha1_hat, "alpha2 =", alpha2_hat, "t1 =", t1_hat, "t2 =", t2_hat, "\n")
  }
  
  LL = sum( log( (alpha1_hat * dexp(data, t1_hat ,log = F)) + (alpha2_hat * dexp(data, t2_hat,log = F)) ) )
  return(c(alpha1_hat,alpha2_hat,t1_hat,t2_hat,LL))
}


dtexp <- function(x, rate=1, lower_trunc, upper_trunc){
  dexp(x, rate, log=T) - logspace.sub(pexp(lower_trunc, rate, lower=F, log=T),pexp(upper_trunc, rate, lower=F, log=T))}


fnx_SP_seg = function(par) -sum(dtexp(l/100, par[1] - offset_vector_segment, lower_trunc=lower_trunc/100,upper_trunc = upper_trunc/100))

# extended pulse on segments 

dtlomax <- function(x, scale, shape, lower_trunc, upper_trunc){
  VGAM::dlomax(x, shape3.q = shape,scale = scale, log=T) -
    logspace.sub(VGAM::plomax(lower_trunc, scale = scale, shape3.q = shape, lower=F, log=T),VGAM::plomax(upper_trunc, scale = scale, shape3.q = shape, lower=F, log=T))}

fnx_EP_seg = function(par){
  tm = par[1] - offset_vector_segment
  td = par[2]
  k = 1/((td/(4*(tm)))^2)
  LL = sum(dtlomax(l/100, scale=((k)/(tm)), shape=((k)+1), lower_trunc=lower_trunc/100,upper_trunc=upper_trunc/100))
}


Get_CI_Dating_fn <- function(est,n_data){
  org <- est
  lwr_approx <- est*(1-(1.96/sqrt(n_data)))
  upr_approx <- est*(1+(1.96/sqrt(n_data)))
  param_CI <- c(org,lwr_approx,upr_approx)
  return(param_CI)
}

LRT <- function(LL_H0, LL_H1,n_param_H0,n_param_H1) {
  ## LL needs to be positive
  lr <- 2 * (LL_H1 -LL_H0)
  p_value <- 1 - pchisq(lr, df = n_param_H1 - n_param_H0)
  
  return(list(LR = lr, p_value = p_value))
}

get_geom_smooth_like_lm_CI <- function(intercept_est,slope_est,data,x,y,p,t_val=NULL){
  n <- nrow(data)
  
  # Calculate predicted values
  data$predicted <- intercept_est + slope_est * x
  
  # Calculate residuals 
  data$residuals <- y - data$predicted
  
  
  # Calculate residual variance
  residual_variance <- sum(data$residuals^2) / (n - p)
  
  # Compute design matrix X
  X <- cbind(1, x)
  
  # Compute inverse of X'X
  XtX_inverse <- solve(t(X) %*% X)
  
  # Compute variance-covariance matrix
  vcov_manual <- residual_variance * XtX_inverse
  
  # Calculate standard errors for each predicted value
  data$residual_sd <- sqrt(rowSums((cbind(1, x) %*% vcov_manual) * cbind(1, x)))
  
  
  # Calculate t-values or take fomr input
  if(is.null(t_val)){
    t_val <- qt(0.975, df = n - p)
  } else {
    t_val <- t_val
  }
  
  
  
  # Calculate margin of error for each predicted value
  data$margin_error <- t_val * data$residual_sd
  
  # Calculate upper and lower bounds of confidence intervals for each predicted value
  data$upper_ci <- data$predicted + data$margin_error
  data$lower_ci <- data$predicted - data$margin_error
  
  return(data)
}

### Neandertal ancestry throughout the genome fn

intersect_with_deserts_fn <- function(desert, est, expand_,end_ = "end",start_ = "start"){
  overlap_T = est %>% filter(start<desert[[end_]] + expand_,
                             end > desert[[start_]]  - expand_,
                             chrom == desert$chrom)
  
  
  overlap_T %>% mutate(desertID = paste0(desert$chrom,":",desert[[start_]],"-",desert[[end_]]),
                       desert_start = desert[[start_]],
                       desert_end = desert[[end_]],
                       desert_start_V = desert$pos_Vernot_2016,
                       desert_end_V = desert$pos_end_Vernot_2016,
                       desert_start_S = desert$pos_Sankararaman_2016,
                       desert_end_S = desert$pos_end_Sankararaman_2016)
  
  
}

### Fig 4 function

fn_plot_frequ_over_time <- function(All_rle_annotated,sample_list,region_chrom,region_start,region_end,col_hex,plot_title=NULL,y_lable="Frequency of Neandertal ancestry",p_size=5,title_size=12,x_axis_size=6){
  xx <- All_rle_annotated %>% 
    filter(chrom == region_chrom,
           pos_end > region_start,
           pos < region_end) %>% 
    filter(SGDP_Superpop %in% c("WEurAs","CeAsSi")) %>% dplyr::select(sample,time_slice) %>% distinct(sample,time_slice, .keep_all = T) %>% mutate(present = T) %>% 
    full_join(.,sample_list,by=c("sample"="sample_name","time_slice"="time_slice")) %>% 
    mutate(present = ifelse(is.na(present),F,present)) %>% 
    mutate(plot_position = ifelse(time_slice == "0","4",ifelse(time_slice == "1","3",ifelse(time_slice == "10000","2","1"))))
  
  
  P <- ggplot() +
    geom_bar(data=xx[xx$time_slice == "0",],aes(x = plot_position,fill = present),position = "fill") +
    geom_point(data=xx[xx$time_slice == "1",],
               aes(x=plot_position,y=seq(from = 0,to = 1,length.out = nrow(xx[xx$time_slice == "1",])),color=present),size = p_size) +
    geom_point(data=xx[xx$time_slice == "10000",],
               aes(x=plot_position,y=seq(from = 0,to = 1,length.out = nrow(xx[xx$time_slice == "10000",])),color=present),size = p_size) +
    geom_point(data=xx[xx$time_slice == "30000",],
               aes(x=plot_position,y=seq(from = 0,to = 1,length.out = nrow(xx[xx$time_slice == "30000",])),color=present),size = p_size) +
    scale_x_discrete(labels=c('45 - 30 ky\nn = 18', '30 - 10 ky\nn = 10', '10000 - 1 ky\nn = 17', 'Present Day\nn = 101')) +
    THEME +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = x_axis_size),
          plot.title = element_text(size = title_size)) +
    xlab("") +
    ylab(y_lable) +
    scale_fill_manual(values=c("#CDC0B0", col_hex),name="Annotation") +
    scale_color_manual(values=c("#CDC0B0", col_hex),name="Annotation") 
  
  if(is.null(plot_title)){
    P + ggtitle(paste0("chr",region_chrom,":",region_start," - ",region_end))
  } else{
    P + ggtitle(plot_title)
  }
  
}    


