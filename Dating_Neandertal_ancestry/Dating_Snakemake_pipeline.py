from collections import Counter
from os import path
import pandas as pd
import numpy as np

configfile: "config/Dating.yml"

####
Folder_name = "Dating_Results"
set_gen_time = 28
###

### Single Dating Sp and EP ###

rule Single_fit_ACov_EP:
    input:
        meta_data = "path_to/Joint_Meta_data_genetic_clusters_using_1240k.csv"
    params:
        gen_time = set_gen_time,
        step_size_td = 1,
        step_size_tm = 1,
        max_td = 300,
        max_tm = 300,
    output:
        output_SP = "{Folder_name}/EarlOoA_single_sample/Single_ACov_EP_SP_fits_samples_{sample_set}_truncation_{truncation}_map_{Map}_min_l_{min_l_ACov}_SP.csv",
        output_EP = "{Folder_name}/EarlOoA_single_sample/Single_ACov_EP_SP_fits_samples_{sample_set}_truncation_{truncation}_map_{Map}_min_l_{min_l_ACov}_EPnlsfit.csv",
        output_EP_grid = "{Folder_name}/EarlOoA_single_sample/Single_ACov_EP_SP_fits_samples_{sample_set}_truncation_{truncation}_map_{Map}_min_l_{min_l_ACov}_EPgridfit.csv", 
        output_summary =  "{Folder_name}/EarlOoA_single_sample/Single_ACov_EP_SP_fits_samples_{sample_set}_truncation_{truncation}_map_{Map}_min_l_{min_l_ACov}_Summary.csv"      
    script:
        "Single_sample_ACov_SP_EP_snake.R"

def merge_all_res_files_SingleSample_fit_ACov_EP_l(wc):
    sample_set = config['set'][wc.set_single_sample_SP_EP_ACov]['sample_set']
    sample_set = ['-'.join(sub_list) for sub_list in sample_set]
    min_l_ACov = config['set'][wc.set_single_sample_SP_EP_ACov]['min_l_ACov']
    weighted = config['set'][wc.set_single_sample_SP_EP_ACov]['weighted']
    Map = config['set'][wc.set_single_sample_SP_EP_ACov]['Map']
    Asc = config['set'][wc.set_single_sample_SP_EP_ACov]['Asc']
    truncation = config['set'][wc.set_single_sample_SP_EP_ACov]['truncation']

    l = expand("{Folder_name}/EarlOoA_single_sample/Single_ACov_EP_SP_fits_samples_{sample_set}_truncation_{truncation}_map_{Map}_min_l_{min_l_ACov}_Summary.csv" ,
        Folder_name=Folder_name,Map=Map,truncation=truncation,sample_set=sample_set,min_l_ACov=min_l_ACov,Asc=Asc,weighted=weighted)
    return l


rule merge_all_res_files_ACov_Single_sample_EP_SP:
    input:
        merge_all_res_files_SingleSample_fit_ACov_EP_l
    output:
        "{Folder_name}/EarlOoA_single_sample/All_ACov_singleSample_SP_EP_summary_{set_single_sample_SP_EP_ACov}.csv"
    shell:
        "head -qn1 {input[0]}  > {output} && "
        "tail -qn+2 {input}  >> {output} "

rule Single_fit_ACov_EP_JKchrom:
    input:
        meta_data = "/home/leonardo_iasi/EMH_Introgression_Project/Introgression_Detection/Joint_Meta_data_genetic_clusters_using_1240k.csv"
    params:
        gen_time = set_gen_time,
        step_size_td = 1,
        step_size_tm = 1,
        max_td = 300,
        max_tm = 300,
    output:
        output_SP = "{Folder_name}/EarlOoA_single_sample/JKchrom/Single_ACov_EP_SP_fits_samples_{sample_set}_truncation_{truncation}_map_{Map}_min_l_{min_l_ACov}_SP_JKchrom_{JKchrom}.csv",
        output_EP = "{Folder_name}/EarlOoA_single_sample/JKchrom/Single_ACov_EP_SP_fits_samples_{sample_set}_truncation_{truncation}_map_{Map}_min_l_{min_l_ACov}_EPnlsfit_JKchrom_{JKchrom}.csv",
        output_EP_grid = "{Folder_name}/EarlOoA_single_sample/JKchrom/Single_ACov_EP_SP_fits_samples_{sample_set}_truncation_{truncation}_map_{Map}_min_l_{min_l_ACov}_EPgridfit_JKchrom_{JKchrom}.csv", 
        output_summary =  "{Folder_name}/EarlOoA_single_sample/JKchrom/Single_ACov_EP_SP_fits_samples_{sample_set}_truncation_{truncation}_map_{Map}_min_l_{min_l_ACov}_Summary_JKchrom_{JKchrom}.csv"    
    script:
        "Single_Sample_ACov_SP_EP_JKchrom_snake.R"

def merge_all_res_files_SingleSample_fit_ACov_EP_JKchrom_l(wc):
    sample_set = config['set'][wc.set_single_sample_SP_EP_ACov_JKchrom]['sample_set']
    sample_set = ['-'.join(sub_list) for sub_list in sample_set]
    min_l_ACov = config['set'][wc.set_single_sample_SP_EP_ACov_JKchrom]['min_l_ACov']
    weighted = config['set'][wc.set_single_sample_SP_EP_ACov_JKchrom]['weighted']
    Map = config['set'][wc.set_single_sample_SP_EP_ACov_JKchrom]['Map']
    Asc = config['set'][wc.set_single_sample_SP_EP_ACov_JKchrom]['Asc']
    truncation = config['set'][wc.set_single_sample_SP_EP_ACov_JKchrom]['truncation']
    n_JKchrom = config['set'][wc.set_single_sample_SP_EP_ACov_JKchrom]['JKchrom']
    JKchrom = range(1,n_JKchrom)

    l = expand("{Folder_name}/EarlOoA_single_sample/JKchrom/Single_ACov_EP_SP_fits_samples_{sample_set}_truncation_{truncation}_map_{Map}_min_l_{min_l_ACov}_Summary_JKchrom_{JKchrom}.csv" ,
        Folder_name=Folder_name,Map=Map,truncation=truncation,sample_set=sample_set,min_l_ACov=min_l_ACov,Asc=Asc,weighted=weighted,JKchrom=JKchrom)
    return l


rule merge_all_res_files_ACov_Singlesample_EP_SP_JKchrom:
    input:
        merge_all_res_files_SingleSample_fit_ACov_EP_JKchrom_l
    output:
        "{Folder_name}/EarlOoA_single_sample/JKchrom/All_ACov_singleSample_SP_EP_summary_{set_single_sample_SP_EP_ACov_JKchrom}.csv"
    shell:
        "head -qn1 {input[0]}  > {output} && "
        "tail -qn+2 {input}  >> {output} "

### Joint Dating ###

""" ACOV """

rule fit_ACov_EP:
    input:
        meta_data = "path_to/Joint_Meta_data_genetic_clusters_using_1240k.csv"
    params:
        gen_time = set_gen_time,
        step_size_td = 1,
        step_size_tm = 1,
        max_td = 2000,
        max_tm = 500,
    output:
        output_LL_matrix = "{Folder_name}/Fit_Data/ACov_EP_fit_{truncation}_map_{Map}_exclude_{exclude_samples}_min_l_{min_l_ACov}_asc_{Asc}_weighted_{weighted}.csv"

    script:
        "ACov_fit_EP_snake.R"

def merge_all_res_files_ACov_EP_l(wc):
    exclude_samples = config['set'][wc.set_ACov_EP]['exclude_samples']
    exclude_samples = ['-'.join(sub_list) for sub_list in exclude_samples]
    min_l_ACov = config['set'][wc.set_ACov_EP]['min_l_ACov']
    weighted = config['set'][wc.set_ACov_EP]['weighted']
    Map = config['set'][wc.set_ACov_EP]['Map']
    Asc = config['set'][wc.set_ACov_EP]['Asc']
    truncation = config['set'][wc.set_ACov_EP]['truncation']

    l = expand("{Folder_name}/Fit_Data/ACov_EP_fit_{truncation}_map_{Map}_exclude_{exclude_samples}_min_l_{min_l_ACov}_asc_{Asc}_weighted_{weighted}.csv",
        Folder_name=Folder_name,Map=Map,truncation=truncation,exclude_samples=exclude_samples,min_l_ACov=min_l_ACov,Asc=Asc,weighted=weighted)
    return l

rule merge_all_res_files_ACov_EP:
    input:
        merge_all_res_files_ACov_EP_l
    output:
        "{Folder_name}/Fit_Data/Results/All_ACov_EP_res_merged_{set_ACov_EP}.csv"
    shell:
        "head -qn1 {input[0]}  > {output} && "
        "tail -qn+2 {input}  >> {output} "

rule fit_ACov_SP:
    input:
        meta_data = "path_to/Joint_Meta_data_genetic_clusters_using_1240k.csv"
    params:
        gen_time = set_gen_time,
    output:
        output_SP = "{Folder_name}/Fit_Data/ACov_SP_fit_{truncation}_map_{Map}_exclude_{exclude_samples}_min_l_{min_l_ACov}_asc_{Asc}_weighted_{weighted}.csv"
    script:
        "ACov_fit_SP_snake.R"

def merge_all_res_files_ACov_SP_l(wc):
    exclude_samples = config['set'][wc.set_ACov_SP]['exclude_samples']
    exclude_samples = ['-'.join(sub_list) for sub_list in exclude_samples]
    min_l_ACov = config['set'][wc.set_ACov_SP]['min_l_ACov']
    weighted = config['set'][wc.set_ACov_SP]['weighted']
    Map = config['set'][wc.set_ACov_SP]['Map']
    Asc = config['set'][wc.set_ACov_SP]['Asc']
    truncation = config['set'][wc.set_ACov_SP]['truncation']

    l = expand("{Folder_name}/Fit_Data/ACov_SP_fit_{truncation}_map_{Map}_exclude_{exclude_samples}_min_l_{min_l_ACov}_asc_{Asc}_weighted_{weighted}.csv",
        Folder_name=Folder_name,Map=Map,truncation=truncation,exclude_samples=exclude_samples,min_l_ACov=min_l_ACov,Asc=Asc,weighted=weighted)
    return l

rule merge_all_res_files_ACov_SP:
    input:
        merge_all_res_files_ACov_SP_l
    output:
        "{Folder_name}/Fit_Data/Results/All_ACov_SP_res_merged_{set_ACov_SP}.csv"
    shell:
        "head -qn1 {input[0]}  > {output} && "
        "tail -qn+2 {input}  >> {output} "


""" JACKKNIFE CHROMS ACOV """

rule fit_EP_ACov_JKchrom:
    input:
        meta_data = "path_to/Joint_Meta_data_genetic_clusters_using_1240k.csv"
    params:
        gen_time = set_gen_time,
        step_size_td = 1,
        step_size_tm = 1,
        max_td = 1000,
        max_tm = 500,
    output:
        output_LL_matrix = "{Folder_name}/Data_JKchrom/ACov_EP_fit_{truncation}_map_{Map}_{exclude_samples}_{min_l_ACov}_{Asc}_{weighted}_JKchrom_{JKchrom}.csv"
    script:
        "ACov_fit_JKchrom_EP_snake.R"


def merge_all_res_files_ACov_EP_JKchrom_l(wc):
    exclude_samples = config['set'][wc.set_ACov_EP_JKchrom]['exclude_samples']
    exclude_samples = ['-'.join(sub_list) for sub_list in exclude_samples]
    min_l_ACov = config['set'][wc.set_ACov_EP_JKchrom]['min_l_ACov']
    weighted = config['set'][wc.set_ACov_EP_JKchrom]['weighted']
    Map = config['set'][wc.set_ACov_EP_JKchrom]['Map']
    Asc = config['set'][wc.set_ACov_EP_JKchrom]['Asc']
    truncation = config['set'][wc.set_ACov_EP_JKchrom]['truncation']
    n_JKchrom = config['set'][wc.set_ACov_EP_JKchrom]['JKchrom']
    JKchrom = range(1,n_JKchrom)

    l = expand( "{Folder_name}/Data_JKchrom/ACov_EP_fit_{truncation}_map_{Map}_{exclude_samples}_{min_l_ACov}_{Asc}_{weighted}_JKchrom_{JKchrom}.csv",
        Folder_name=Folder_name,Map=Map,truncation=truncation,exclude_samples=exclude_samples,min_l_ACov=min_l_ACov,Asc=Asc,weighted=weighted,JKchrom=JKchrom)
    return l


rule merge_all_res_files_ACov_EP_JKchrom:
    input:
        merge_all_res_files_ACov_EP_JKchrom_l
    output:
        "{Folder_name}/Data_JKchrom/Results/All_ACov_EP_res_merged_{set_ACov_EP_JKchrom}.csv"
    shell:
        "head -qn1 {input[0]}  > {output} && "
        "tail -qn+2 {input}  >> {output} "


""" Segments """

rule fit_Seg_SP:
    input:
        meta_data = "path_to/Joint_Meta_data_genetic_clusters_using_1240k.csv"
    params:
        gen_time = set_gen_time
    output:
        output_SP = "{Folder_name}/Fit_Data/Seg_SP_fit_map_{Map}_exclude_{exclude_samples}_min_l_{min_l_Seg}_asc_{Asc}.csv"

    script:
        "Seg_fit_SP_snake.R"

def merge_all_res_files_Seg_SP_l(wc):
    exclude_samples = config['set'][wc.set_Seg_SP]['exclude_samples']
    exclude_samples = ['-'.join(sub_list) for sub_list in exclude_samples]
    min_l_Seg = config['set'][wc.set_Seg_SP]['min_l_Seg']
    Map = config['set'][wc.set_Seg_SP]['Map']
    Asc = config['set'][wc.set_Seg_SP]['Asc']


    l = expand("{Folder_name}/Fit_Data/Seg_SP_fit_map_{Map}_exclude_{exclude_samples}_min_l_{min_l_Seg}_asc_{Asc}.csv",
        Folder_name=Folder_name,Map=Map,exclude_samples=exclude_samples,min_l_Seg=min_l_Seg,Asc=Asc)
    return l

rule merge_all_res_files_Seg_SP:
    input:
        merge_all_res_files_Seg_SP_l
    output:
        "{Folder_name}/Fit_Data/Results/All_Seg_SP_res_merged_{set_Seg_SP}.csv"
    shell:
        "head -qn1 {input[0]}  > {output} && "
        "tail -qn+2 {input}  >> {output} "


def exec_admixfrog(wc):
  
    set_ACov_SP = "ACov_Run_SP"
    set_Seg_SP = "Seg_Run_SP"

    set_single_sample_SP_EP_ACov = "Single_Sample_EP_SP"
    set_single_sample_SP_EP_ACov_JKchrom = "Single_Sample_EP_SP_JKchrom"

    set_ACov_EP = "ACov_Run_EP"
    set_ACov_EP_JKchrom = "ACov_Run_EP_JKchrom"



    ## SP fits

    file_1 = expand("{Folder_name}/Fit_Data/Results/All_ACov_SP_res_merged_{set_ACov}.csv",Folder_name=Folder_name,set_ACov=set_ACov_SP)

    file_2 = expand("{Folder_name}/Fit_Data/Results/All_Seg_SP_res_merged_{set_Seg}.csv",Folder_name=Folder_name,set_Seg=set_Seg_SP)

    ## EP single sample

    file_3 = expand("{Folder_name}/EarlOoA_single_sample/All_ACov_singleSample_SP_EP_summary_{set_single_sample_SP_EP_ACov}.csv",Folder_name=Folder_name,set_single_sample_SP_EP_ACov=set_single_sample_SP_EP_ACov)
    file_3_1 = expand("{Folder_name}/EarlOoA_single_sample/JKchrom/All_ACov_singleSample_SP_EP_summary_{set_single_sample_SP_EP_ACov_JKchrom}.csv",Folder_name=Folder_name,set_single_sample_SP_EP_ACov_JKchrom=set_single_sample_SP_EP_ACov_JKchrom)

    ## EP joint

    file_4 = expand("{Folder_name}/Fit_Data/Results/All_ACov_EP_res_merged_{set_ACov}.csv",Folder_name=Folder_name,set_ACov=set_ACov_EP)
    file_4_1 = expand("{Folder_name}/Data_JKchrom/Results/All_ACov_EP_res_merged_{set_ACov_JKchrom}.csv",Folder_name=Folder_name,set_ACov_JKchrom=set_ACov_EP_JKchrom)


    files = file_1 + file_2 + file_3 + file_3 + file_3_1 + file_4 + file_4_1
    return files

rule run_series_cfg:
    input:
        exec_admixfrog
    output: touch('controller/Fit')