### functions
source("~/EMH_Introgression_Project/Introgression_Detection/scripts/Analysis_R_functions.R")

sample_two_pulse_ACov_sim <- function(min_dist,max_dist,binwidth,A1,A2,tm1,tm2,c=0,noise=0.001,offset,trunc_,round_dec=0.01){
  dist= seq(min_dist,max_dist,binwidth)
  wcov = c + A1*exp(-dist*(tm1-offset)) + A2*exp(-dist*(tm2-offset))
  wcov_n = wcov + rnorm(length(wcov),mean = 0,sd = noise)
  data.x = data.frame(dist=dist,wcov=wcov_n) 
  if(trunc_ == "estimated_truncation"){
    max_l  = data.x$dist[which(data.x$wcov <= 1e-6)[1]]
    max_l = round_any(max_l,round_dec, f = ceiling) 
  } else {
    trunc_ = as.numeric(trunc_)
    max_l = trunc_
  }
  data.x = data.x %>% filter(dist <= max_l)
  if(weighted == "yes"){
    data.x <- data.x %>% mutate(weight = (max_dist / binwidth) /length(dist))
  } else {
    data.x <- data.x %>% mutate(weight = 1)
  }
  return(data.x)
}

sample_two_pulse_Seg_sim <- function(n_seg,tm1,tm2,offset,noise){
  data.x = rexp(n_seg,rate = (tm1 - offset)) + rexp(n_seg,rate = (tm2 - offset))
  data.x_n = data.x + rnorm(length(data.x),mean = 0,sd = noise/100)
  return(data.x_n)
}


### input Snakemake

Meta = snakemake@input$meta_data

tm1 = as.numeric(snakemake@wildcards$tm1)
tm2 = as.numeric(snakemake@wildcards$tm2)
trunc_ = snakemake@wildcards$truncation
min_cut_ACov = as.numeric(snakemake@wildcards$min_ACov)
min_len_Seg = as.numeric(snakemake@wildcards$min_Seg)
noise = as.numeric(snakemake@wildcards$noise)
n_rep = snakemake@wildcards$rep

output = snakemake@output$output

if(trunc_ == "estimated_truncation"){
  weighted = "yes"
} else {
  weighted = "no"
}

## ACov
binwidth = 0.00005
min_dist = 0
max_dist = 0.1
max_M = 0.1
true_A1_range = c(0.02,0.04)
true_A2_range = c(0.01,0.02)
true_C_range = c(0,0.002)



## segments
min_len = 0
n_seg = 150


# Meta
Joint_Meta <- read.csv(Meta) %>% mutate(PCA_pops = ifelse(ML_BP_Mean == 0 ,pop,sample_name),Cov = ifelse(ML_BP_Mean == 0 ,ave_cov_reported,ave_cov_array),DataType = ifelse(Data_source == "Shotgun" | Data_source == "Genotypes","Shotgun","Captured" )) %>% filter(SGDP_Superpop != "Africa" )

gen_time = 29
min_age_BP = 20000

n_inds = length(Joint_Meta$sample_name[Joint_Meta$ML_BP_Mean >= min_age_BP])
offset = Joint_Meta$ML_BP_Mean[Joint_Meta$ML_BP_Mean >= min_age_BP] / gen_time

# define th grid
td_values = c(1,seq(from = 10, to = 2000, by = 10))
tm_values = c(round(max(offset)) + 1, round(max(offset)) + seq(from = 10, to = 500, by = 10))

true_A1 = runif(n_inds,true_A1_range[1],true_A1_range[2])
true_A2 = runif(n_inds,true_A2_range[1],true_A2_range[2])
true_C = runif(n_inds,true_C_range[1],true_C_range[1])

sample.list.TP <- list()
for(i in 1:length(offset)){
  sample.list.TP[[i]] <- sample_two_pulse_ACov_sim(min_dist=min_dist,max_dist=max_dist,binwidth=binwidth
                                                   ,A1 =true_A1[i],A2 =true_A2[i],tm1 = tm1,tm2 = tm2,c=true_C[i],noise = noise,offset  = offset[i],
                                                   trunc_  = trunc_)
}


### simulated Segment data ###

sample.list.TP_seg_l <- list()
for(i in 1:length(offset)){
  x <- sample_two_pulse_Seg_sim(n_seg = n_seg,tm1 = tm1,tm2 = tm2,offset = offset[i],noise=noise) 
  x <- x[x >= min_len/100]
  sample.list.TP_seg_l[[i]] <- data.frame(seg=x,ts=rep(offset[i],length(x )))
}

## ACov
sample.list.TP_trunc <- list()
for(i in 1:length(sample.list.TP)){
  sample.list.TP_trunc[[i]] <- sample.list.TP[[i]] %>% filter(dist >= min_cut_ACov)
}



data.in = list(
  dist=lapply(sample.list.TP_trunc, `[[`, 'dist'),
  ACov_res_list=sample.list.TP_trunc,
  y=unlist(lapply(sample.list.TP_trunc, `[[`, 'wcov')),
  weight = unlist(lapply(sample.list.TP_trunc, `[[`, 'weight')),
  offset = offset,
  A_init = rep(0.1,length(offset)),
  tm_init = max(offset)+50,
  C_init = rep(0,length(offset))
)


res_SP = Joint_Simple_Pulse_ACov_weighted_fn(data.in = data.in,minFac = 1e-6,maxIter = 10,JointInterceptsFit = F)

tm_CI_SP_ACov <- Get_CI_Dating_fn(est = res_SP$tm,n_data = length(data.in$y))
RSS_SP_ACov <- res_SP$RSS

RSS_SP_ACov <- res_SP$RSS
n_DP = length(data.in$y)
LL_SP_ACov_fit = - n_DP/2 * log(RSS_SP_ACov) - n_DP/2 * log(2*pi / n_DP) - n_DP/2

Sp_ACov_res <- data.frame(tm = tm_CI_SP_ACov[1],tm_lower = tm_CI_SP_ACov[2],tm_upper = tm_CI_SP_ACov[3],LL=LL_SP_ACov_fit,
                          method="ACov_fit",model="SimplePulse",simulation="TwoPulses",min_l = min_cut_ACov,truncation = trunc_,weighted=weighted,
                          true_tm1= tm1,true_tm2= tm2,true_td = 1,n_rep = n_rep,noise=noise)

## segments

sample.list.TP_seg_l_trunc <- list()
for(i in 1:length(sample.list.TP_seg_l)){
  sample.list.TP_seg_l_trunc[[i]] <- sample.list.TP_seg_l[[i]] %>% filter(seg >= min_len_Seg/100)
}

l = unlist(lapply(sample.list.TP_seg_l_trunc, `[[`, 'seg'))*100
offset_vector_segment = unlist(lapply(sample.list.TP_seg_l_trunc, `[[`, 'ts'))
lower_trunc = min_len_Seg
upper_trunc = max(l)
tm_lower = max(offset_vector_segment) + 1
tm_upper = max(offset_vector_segment) + 1000

res_seg_SP = optim(c(tm_lower + 10), fnx_SP_seg, method="L-BFGS-B",lower = c(tm_lower),upper = c(tm_upper))
tm_SP_seg_fit <- round(Get_CI_Dating_fn(est = res_seg_SP$par,n_data = length(l)))
LL_SP_seg_fit = res_seg_SP$value

Sp_seg_res <- data.frame(tm = tm_SP_seg_fit[1],tm_lower = tm_SP_seg_fit[2],tm_upper = tm_SP_seg_fit[3],LL=LL_SP_seg_fit,
                         method="Segments",model="SimplePulse",simulation="TwoPulses",min_l = min_len,truncation = NA,weighted=NA,
                         true_tm1= tm1,true_tm2= tm2,true_td = 1,n_rep = n_rep,noise=noise)



## EP fit


### ACon fit
res_matrix_ACov_fit <- matrix(NA,nrow = length(td_values),ncol =  length(tm_values),byrow = T)
res_amplitude.SP <- matrix(NA,nrow = length(td_values) * length(tm_values),ncol =  length(offset),byrow = T)
res_C.SP <- matrix(NA,nrow = length(td_values) * length(tm_values),ncol =  length(offset),byrow = T)


for(i in 1:length(td_values)){
  for(j in 1:length(tm_values)){
    res_fit <- GriD_EP_ACov_amplitude_x_inter_fit_weighted(data.in = data.in,td_fixed = td_values[i],tm_fixed = tm_values[j],A_init = data.in$A_init,C_init = data.in$C_init)
    res_matrix_ACov_fit[i,j] <- res_fit$RSS
    res_amplitude.SP[((i+j)/2),] <- res_fit$A
    res_C.SP[((i+j)/2),] <- res_fit$C
  }
}

colnames(res_matrix_ACov_fit) <- tm_values
rownames(res_matrix_ACov_fit) <- td_values

n_DP = length(data.in$y)
res_matrix_ACov_fit_LL <- matrix(NA,nrow = length(td_values),ncol =  length(tm_values),byrow = T)

for(i in 1:length(td_values)){
  for(j in 1:length(tm_values)){
    res_matrix_ACov_fit_LL[i,j] <-  - n_DP/2 * log(res_matrix_ACov_fit[i,j]) - n_DP/2 * log(2*pi / n_DP) - n_DP/2
    
  }
}

colnames(res_matrix_ACov_fit_LL) <- tm_values
rownames(res_matrix_ACov_fit_LL) <- td_values

# ACov fit
res_matrix_ACov_fit_LL <- melt(res_matrix_ACov_fit_LL)
res_matrix_ACov_fit_LL$value0 <- res_matrix_ACov_fit_LL$value - max(res_matrix_ACov_fit_LL$value)

est_td_EP_ACov_fit_EP=res_matrix_ACov_fit_LL$Var1[res_matrix_ACov_fit_LL$value0 == max(res_matrix_ACov_fit_LL$value0)]
est_tm_EP_ACov_fit_EP=res_matrix_ACov_fit_LL$Var2[res_matrix_ACov_fit_LL$value0 == max(res_matrix_ACov_fit_LL$value0)]

res_matrix_ACov_fit_LL$value0 <- ifelse(res_matrix_ACov_fit_LL$value0 < -100,-101,res_matrix_ACov_fit_LL$value0)
res_matrix_ACov_fit_LL$method <- "ACov_fit"
res_matrix_ACov_fit_LL$model <- "ExtendedPulse"
res_matrix_ACov_fit_LL$simulation <- "TwoPulses"
res_matrix_ACov_fit_LL$est_tm <- est_tm_EP_ACov_fit_EP
res_matrix_ACov_fit_LL$est_td <- est_td_EP_ACov_fit_EP

res_matrix_ACov_fit_LL$true_tm1 <-  tm1
res_matrix_ACov_fit_LL$true_tm2 <-  tm2
res_matrix_ACov_fit_LL$true_td <-  1
res_matrix_ACov_fit_LL$n_rep <- n_rep

res_matrix_ACov_fit_LL$min_l <- min_len
res_matrix_ACov_fit_LL$weighted <- weighted
res_matrix_ACov_fit_LL$truncation <- trunc_
res_matrix_ACov_fit_LL$noise <- noise

# segments

res_matrix_Segfit <- matrix(NA,nrow = length(td_values),ncol =  length(tm_values),byrow = T)

for(i in 1:length(td_values)){
  for(j in 1:length(tm_values)){
    res_matrix_Segfit[i,j] <- fnx_EP_seg(c(tm_values[j],td_values[i]))
    
  }
}

colnames(res_matrix_Segfit) <- tm_values
rownames(res_matrix_Segfit) <- td_values

res_matrix_Segfit <- reshape2::melt(res_matrix_Segfit)
res_matrix_Segfit$value0 <-  res_matrix_Segfit$value - max(res_matrix_Segfit$value)


est_td_EP_seg_fit_EP=res_matrix_Segfit$Var1[res_matrix_Segfit$value0 == max(res_matrix_Segfit$value0)]
est_tm_EP_seg_fit_EP=res_matrix_Segfit$Var2[res_matrix_Segfit$value0 == max(res_matrix_Segfit$value0)]

res_matrix_Segfit$value0 <- ifelse(res_matrix_Segfit$value0 < -100,-101,res_matrix_Segfit$value0)
res_matrix_Segfit$method <- "Segments"
res_matrix_Segfit$model <- "ExtendedPulse"
res_matrix_Segfit$est_tm <- est_tm_EP_seg_fit_EP
res_matrix_Segfit$est_td <- est_td_EP_seg_fit_EP

res_matrix_Segfit$true_tm1 <-  tm1
res_matrix_Segfit$true_tm2 <-  tm2
res_matrix_Segfit$true_td <-  1
res_matrix_Segfit$n_rep <- n_rep

res_matrix_Segfit$min_l <- min_len
res_matrix_Segfit$simulation <- "TwoPulses"
res_matrix_Segfit$weighted <- NA
res_matrix_Segfit$truncation <- NA
res_matrix_Segfit$noise <- noise



save(Sp_ACov_res,Sp_seg_res,res_matrix_ACov_fit_LL,res_matrix_Segfit,file = output)
