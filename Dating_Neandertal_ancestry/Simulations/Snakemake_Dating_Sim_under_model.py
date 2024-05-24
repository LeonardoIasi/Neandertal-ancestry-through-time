from collections import Counter
from os import path
import pandas as pd
import numpy as np

configfile: "Dating_Sim.yaml"

####
path_to_Dryad_folder = "my_path"
Folder_name = "Dating_Results_Sim"
####


rule sim_SP:
    input:
        meta_data = "{path_to_Dryad_folder}/Meta_Data_individuals.csv"
    output:
        output= "{Folder_name}/SP_sim_{truncation}_tm_{tm}_minACov_{min_ACov}_minSeg_{min_Seg}_noise_{noise}_rep_{rep}.RData",
    script:
        "SP_sim_snake.R"

rule sim_TP:
    input:
        meta_data = "{path_to_Dryad_folder}/Meta_Data_individuals.csv",
    output:
        output= "{Folder_name}/TP_sim_{truncation}_tm1_{tm1}_tm2_{tm2}_minACov_{min_ACov}_minSeg_{min_Seg}_noise_{noise}_rep_{rep}.RData",
    script:
        "TP_sim_snake.R"


rule sim_EP:
    input:
        meta_data = "{path_to_Dryad_folder}/Meta_Data_individuals.csv",
    output:
        output= "{Folder_name}/EP_sim_{truncation}_tm_td_comb_{tm_td_comb}_minACov_{min_ACov}_minSeg_{min_Seg}_noise_{noise}_rep_{rep}.RData",
    script:
        "EP_sim_snake.R"



def exec_admixfrog(wc):
    sim_Set = "SP_Sim"
    min_l_ACov = config['set'][sim_Set]['min_l_ACov']
    min_l_Seg = config['set'][sim_Set]['min_l_Seg']
    tm = config['set'][sim_Set]['tm']
    td = config['set'][sim_Set]['td']
    noise = config['set'][sim_Set]['noise']
    truncation = config['set'][sim_Set]['truncation']
    n_rep = config['set'][sim_Set]['n_rep']
    rep = range(0,n_rep)

    file_1 = expand("{Folder_name}/SP_sim_{truncation}_tm_{tm}_minACov_{min_ACov}_minSeg_{min_Seg}_noise_{noise}_rep_{rep}.RData",
        Folder_name = Folder_name,tm=tm,truncation=truncation,min_Seg=min_l_Seg,min_ACov=min_l_ACov,noise=noise,rep=rep)

    sim_Set = "TP_Sim"
    tm1 = config['set'][sim_Set]['tm1']
    tm2= config['set'][sim_Set]['tm2']

    file_2 = expand("{Folder_name}/TP_sim_{truncation}_tm1_{tm1}_tm2_{tm2}_minACov_{min_ACov}_minSeg_{min_Seg}_noise_{noise}_rep_{rep}.RData",
        Folder_name = Folder_name,tm1=tm1,tm2=tm2,truncation=truncation,min_Seg=min_l_Seg,min_ACov=min_l_ACov,noise=noise,rep=rep)
    
    sim_Set = "EP_Sim"
    tm_td_comb = config['set'][sim_Set]['tm_td_comb']
    tm_td_comb = ['{}_{}'.format(i[0],i[1]) for i in tm_td_comb]

    file_3 = expand("{Folder_name}/EP_sim_{truncation}_tm_td_comb_{tm_td_comb}_minACov_{min_ACov}_minSeg_{min_Seg}_noise_{noise}_rep_{rep}.RData",
        Folder_name = Folder_name,tm_td_comb=tm_td_comb,truncation=truncation,min_Seg=min_l_Seg,min_ACov=min_l_ACov,noise=noise,rep=rep)

    files = file_1 + file_2 + file_3
    return files

rule run_series_cfg:
    input:
        exec_admixfrog
    output: touch('controller/Sim')