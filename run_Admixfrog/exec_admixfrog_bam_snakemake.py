from collections import Counter
from os import path
import pandas as pd
import numpy as np

configfile: "config/data_example.yml"
configfile: "config/panel_example.yml"

include: "Admixfrog_snakemake_for_bams_F4_calculation.snake"
include: "Admixfrog_snakemake_for_bams_Overview_Stats.snake"
include: "Phylogenetic_fragment_analysis/Phylogenetic_fragment_analysis_fromBAM_geneticMap.snake"
include: "Admixfrog_snakemake_for_bams.snake"

CHROMS = [str(i+1) for i in range(22)] + ["X"]

Folder_name = "Result_Folder"


""" EXECUTE SNAKEMAKE """
def exec_admixfrog(wc):
    Admixfrog_parameter = "error2"
    BinSize = 5000
    States = "-".join(["AFR", "NEA","DEN"])
    Sample = config['panels'][wc.panel]
    Snpset = wc.snpset
    penalty = [0.25,0.4]
    length_bin = 15
    Ref_match_data = ["gtLH"]
    called_Ancestry = ["NEA","DEN"]
    bin_called_min_len = [0.05,0.2]
    bin_called_min_len_pos = [0]
    min_cov = [1]
    min_n_all_SNP = [0]
    min_n_SNP = [0]
    type = ['state','homo']
    DiagSites = ["OnlyVariable","MatchingAscertainment"]
    map = ["AA_Map","Shared_Map"]


    file_0 =  expand("{Folder_name}/samples/length_bin_{LB}/{D}_{E}.in.xz",
    Folder_name=Folder_name,LB=length_bin,D=Sample,E=Snpset)

    file_1 = expand("{Folder_name}/rle/{A}/length_bin_{LB}/{B}/{C}/{M}/{D}_{E}.rle{F}.xz.anno_all_types",
    Folder_name=Folder_name,A=Admixfrog_parameter,LB=length_bin,B=BinSize,C=States,D=Sample,E=Snpset,F=penalty,M=map)

    file_2 = expand("{Folder_name}/stats/{A}/length_bin_{LB}/{M}/{C}_{P}_{B}_{E}_{MC}.proportion_callable_genome.txt",
    Folder_name=Folder_name,A=Admixfrog_parameter,LB=length_bin,B=BinSize,C=States,E=Snpset,P=wc.panel,M=map,MC=min_cov)

    file_3 = expand("{Folder_name}/stats/{A}/{C}/{M}/QC_{P}_length_bin_{LB}_{B}_{E}.wide_updated.csv",
    Folder_name=Folder_name,A=Admixfrog_parameter,LB=length_bin,B=BinSize,C=States,E=Snpset,M=map,P=wc.panel)

    file_4 = expand("{Folder_name}/stats/{A}/length_bin_{LB}/{B}/{C}/{M}/All_bins_{CA}_{P}_{E}.bin_called{F}_min_len{Min}_min_len_pos{MP}_min_n_all_SNP{allSNP}_min_n_SNP{SNP}.xz",
    Folder_name=Folder_name,A=Admixfrog_parameter,LB=length_bin,B=BinSize,C=States,P=wc.panel,E=Snpset,F=penalty,CA=called_Ancestry,Min=bin_called_min_len,MP=bin_called_min_len_pos,allSNP=min_n_all_SNP,SNP=min_n_SNP,M=map)

    file_5 = expand("{Folder_name}/Diagnostic_sites/{A}/length_bin_{LB}/{B}/{C}/{M}/{DiagSites}/{D}_{E}.rle{F}.match.xz",
                    Folder_name=Folder_name,A=Admixfrog_parameter,LB=length_bin,B=BinSize,C=States,D=Sample,E=Snpset,F=penalty,DiagSites=DiagSites,M=map)
    
    file_6 = expand("{Folder_name}/stats/{A}/length_bin_{LB}/{B}/{C}/{M}/All_bins_raw_{CA}_{P}_{E}.bin_raw.xz",
    Folder_name=Folder_name,A=Admixfrog_parameter,LB=length_bin,B=BinSize,C=States,P=wc.panel,E=Snpset,CA=["NEA","DEN","AFR"],M=map)



    files = file_0 + file_1 + file_2 + file_3 + file_4 + file_5 + file_6
    return files

rule run_series_cfg:
    input:
        exec_admixfrog
    output: touch('../controller/{panel}_{snpset}')