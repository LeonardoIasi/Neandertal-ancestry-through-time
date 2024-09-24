from collections import Counter
from os import path
import pandas as pd
import numpy as np

configfile: "../config/data_EMH_Introgression_Project.yaml"
configfile: "../config/panels_EMH_Introgression_Project.yaml"


CHROMS = [str(i+1) for i in range(22)] + ["X"]
Folder_name = '/mnt/diversity/leonardo_iasi/EMH_Introgression_Project/NewTechnicalValidations_newREF_admixfrog0.7'
Result_Folder_name = '/home/leonardo_iasi/EMH_Introgression_Project/Introgression_Detection/NewTechnicalValidations/Checking_Ancestry_info_newREF_admixfrog0.7'
#TruePenalty = [0.2,0.25,0.3,0.4]

rule admixfrog_sample_input:
    input :
        ref = "/home/leonardo_iasi/EMH_Introgression_Project/Introgression_Detection/admixfrogREF/ref/published/EMHPanel/ref_{snpset}_hs37mMask35to99.csv.xz",
        bam = "{Folder_name}/admixfrog_bams/{snpset}/{sample}.bam",
        bai = "{Folder_name}/admixfrog_bams/{snpset}/{sample}.bam.bai",

    params:
        cutoff = 3,
    output:
        csv = "{Folder_name}/samples/length_bin_{length_bin}/{sample}_{snpset}.in.xz"
    shell:
        "admixfrog-bam --ref {input.ref} --bamfile {input.bam} --force-bam"
        " --deam-cutoff {params.cutoff} "
        " --out {output.csv} "
        " --length-bin-size {wildcards.length_bin} "

def _inputfiles_panel(wc):
    l = expand("{Folder_name}/samples/length_bin_{length_bin}/{sample}_{{snpset}}.in.xz",
        Folder_name=Folder_name,sample = wc.sample,length_bin=wc.length_bin)
    return l


rule run_admixfrog_error2:
    input:
        infile = _inputfiles_panel,
        ref = "/home/leonardo_iasi/EMH_Introgression_Project/Introgression_Detection/admixfrogREF/ref/published/EMHPanel/Shared_Map/ref_{snpset}_hs37mMask35to99.csv.xz",
    params:
        cont_id = "AFR",
        ll_tol = 1e-2,
        freq_f = 3,
        freq_c = 3,
        error = 1e-2,
        ancestral = "PAN",
        run_penalty = 0.1,
        max_iter = 250,
        n_post_rep = 200,
    version: ".3"
    output:
        res_bin = "{Folder_name}/admixfrog/error2/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.bin.xz",
        res_snp = "{Folder_name}/admixfrog/error2/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.snp.xz",
        res_cont= "{Folder_name}/admixfrog/error2/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.cont.xz",
        res_par = "{Folder_name}/admixfrog/error2/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.pars.yaml",
        res_rle = "{Folder_name}/admixfrog/error2/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.rle.xz",
        res_res = "{Folder_name}/admixfrog/error2/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.res.xz",
        res2_res = "{Folder_name}/admixfrog/error2/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.res2.xz",
    run:
        outname = path.splitext(path.splitext(output.res_bin)[0])[0]
        states = wildcards.states.split("-")
        downsample = wildcards.downsampled
        genetic_map = str(wildcards.map)
        s = "admixfrog --infile {input.infile} --ref {input.ref} -o {outname} "
        s += " --states {states} "
        s += " --cont-id {params.cont_id} "
        s += " --ll-tol {params.ll_tol} "
        s += " --bin-size {wildcards.bin_size} "
        s += " --est-F --est-tau --freq-F {params.freq_f} "
        s += " --freq-contamination {params.freq_c} "
        s += " --e0 {params.error} "
        s += " --est-error "
        s += " --ancestral {params.ancestral} "
        s += " --run-penalty {params.run_penalty}"
        s += " --max-iter {params.max_iter}"
        s += " --n-post-replicates {params.n_post_rep}"
        s += " --filter-pos 50 --filter-map 0.000"
        s += " --map-column {genetic_map}"
        s += " --init-guess AFR"
        s += " --downsample {downsample}"
        shell(s)

rule run_admixfrog_error2UnpolarizedSNPsRemoved:
    input:
        infile = _inputfiles_panel,
        ref = "/home/leonardo_iasi/EMH_Introgression_Project/Introgression_Detection/admixfrogREF/ref/published/EMHPanel/Shared_Map/ref_{snpset}_hs37mMask35to99.csv.xz",
    params:
        cont_id = "AFR",
        ll_tol = 1e-2,
        freq_f = 3,
        freq_c = 3,
        error = 1e-2,
        ancestral = "PAN",
        run_penalty = 0.1,
        max_iter = 250,
        n_post_rep = 200,
    version: ".3"
    output:
        res_bin = "{Folder_name}/admixfrog/error2UnpolarizedSNPsRemoved/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.bin.xz",
        res_snp = "{Folder_name}/admixfrog/error2UnpolarizedSNPsRemoved/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.snp.xz",
        res_cont= "{Folder_name}/admixfrog/error2UnpolarizedSNPsRemoved/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.cont.xz",
        res_par = "{Folder_name}/admixfrog/error2UnpolarizedSNPsRemoved/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.pars.yaml",
        res_rle = "{Folder_name}/admixfrog/error2UnpolarizedSNPsRemoved/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.rle.xz",
        res_res = "{Folder_name}/admixfrog/error2UnpolarizedSNPsRemoved/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.res.xz",
        res2_res = "{Folder_name}/admixfrog/error2UnpolarizedSNPsRemoved/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.res2.xz",
    run:
        outname = path.splitext(path.splitext(output.res_bin)[0])[0]
        states = wildcards.states.split("-")
        downsample = wildcards.downsampled
        genetic_map = str(wildcards.map)
        s = "admixfrog --infile {input.infile} --ref {input.ref} -o {outname} "
        s += " --states {states} "
        s += " --cont-id {params.cont_id} "
        s += " --ll-tol {params.ll_tol} "
        s += " --bin-size {wildcards.bin_size} "
        s += " --est-F --est-tau --freq-F {params.freq_f} "
        s += " --freq-contamination {params.freq_c} "
        s += " --e0 {params.error} "
        s += " --est-error "
        s += " --ancestral {params.ancestral} "
        s += " --run-penalty {params.run_penalty}"
        s += " --max-iter {params.max_iter}"
        s += " --n-post-replicates {params.n_post_rep}"
        s += " --filter-pos 50 --filter-map 0.000"
        s += " --map-column {genetic_map}"
        s += " --filter-ancestral"
        s += " --init-guess AFR"
        s += " --downsample {downsample}"
        shell(s)

rule run_admixfrog_error2NoAncestralInfo:
    input:
        infile = _inputfiles_panel,
        ref = "/home/leonardo_iasi/EMH_Introgression_Project/Introgression_Detection/admixfrogREF/ref/published/EMHPanel/Shared_Map/ref_{snpset}_hs37mMask35to99.csv.xz",
    params:
        cont_id = "AFR",
        ll_tol = 1e-2,
        freq_f = 3,
        freq_c = 3,
        error = 1e-2,
        run_penalty = 0.1,
        max_iter = 250,
        n_post_rep = 200,
    version: ".3"
    output:
        res_bin = "{Folder_name}/admixfrog/error2NoAncestralInfo/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.bin.xz",
        res_snp = "{Folder_name}/admixfrog/error2NoAncestralInfo/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.snp.xz",
        res_cont= "{Folder_name}/admixfrog/error2NoAncestralInfo/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.cont.xz",
        res_par = "{Folder_name}/admixfrog/error2NoAncestralInfo/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.pars.yaml",
        res_rle = "{Folder_name}/admixfrog/error2NoAncestralInfo/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.rle.xz",
        res_res = "{Folder_name}/admixfrog/error2NoAncestralInfo/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.res.xz",
        res2_res = "{Folder_name}/admixfrog/error2NoAncestralInfo/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.res2.xz",
    run:
        outname = path.splitext(path.splitext(output.res_bin)[0])[0]
        states = wildcards.states.split("-")
        downsample = wildcards.downsampled
        genetic_map = str(wildcards.map)
        s = "admixfrog --infile {input.infile} --ref {input.ref} -o {outname} "
        s += " --states {states} "
        s += " --cont-id {params.cont_id} "
        s += " --ll-tol {params.ll_tol} "
        s += " --bin-size {wildcards.bin_size} "
        s += " --est-F --est-tau --freq-F {params.freq_f} "
        s += " --freq-contamination {params.freq_c} "
        s += " --e0 {params.error} "
        s += " --est-error "
        s += " --run-penalty {params.run_penalty}"
        s += " --max-iter {params.max_iter}"
        s += " --n-post-replicates {params.n_post_rep}"
        s += " --filter-pos 50 --filter-map 0.000"
        s += " --map-column {genetic_map}"
        s += " --init-guess AFR"
        s += " --downsample {downsample}"
        shell(s)

rule run_admixfrog_gtmode:
    input:
        infile = _inputfiles_panel,
        ref = "/home/leonardo_iasi/EMH_Introgression_Project/Introgression_Detection/admixfrogREF/ref/published/EMHPanel/Shared_Map/ref_{snpset}_hs37mMask35to99.csv.xz",
    params:
        cont_id = "AFR",
        ll_tol = 1e-2,
        freq_f = 3,
        freq_c = 3,
        error = 1e-2,
        ancestral = "PAN",
        run_penalty = 0.1,
        max_iter = 250,
        n_post_rep = 200,
    output:
        res_bin = "{Folder_name}/admixfrog/gtmode/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.bin.xz",
        res_snp = "{Folder_name}/admixfrog/gtmode/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.snp.xz",
        res_cont= "{Folder_name}/admixfrog/gtmode/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.cont.xz",
        res_par = "{Folder_name}/admixfrog/gtmode/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.pars.yaml",
        res_rle = "{Folder_name}/admixfrog/gtmode/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.rle.xz",
        res_res = "{Folder_name}/admixfrog/gtmode/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.res.xz",
        res2_res = "{Folder_name}/admixfrog/gtmode/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.res2.xz",
    run:
        outname = path.splitext(path.splitext(output.res_bin)[0])[0]
        states = wildcards.states.split("-")
        downsample = wildcards.downsampled
        genetic_map = str(wildcards.map)
        s = "admixfrog --infile {input.infile} --ref {input.ref} -o {outname} "
        s += " --filter-pos 50 --filter-map 0"
        s += " --female --gt-mode "
        s += " --states {states} "
        s += " --ll-tol {params.ll_tol} "
        s += " --bin-size {wildcards.bin_size} "
        s += " --est-F --est-tau --freq-F {params.freq_f} "
        s += " --dont-est-contamination "
        s += " --c0 0"
        s += " --map-column {genetic_map}"
        s += " --e0 {params.error} "
        s += " --ancestral {params.ancestral} "
        s += " --run-penalty {params.run_penalty}"
        s += " --max-iter {params.max_iter}"
        s += " --n-post-replicates {params.n_post_rep}"
        s += " --downsample {downsample}"
        shell(s)

rule call_runs:
    input:
        rle = "{Folder_name}/admixfrog/{admf_params}/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}_downsampled.{downsampled}.bin.xz",
    params:
        run_penalty = 0.25,
    output:
       rle = "{Folder_name}/rle/{admf_params}/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}.rle{penalty}_downsampled.{downsampled}.xz",
    shell:
        "admixfrog-rle --out {output.rle} --in {input.rle} "
        "--run-penalty {wildcards.penalty} "

""" BAM PREPARATION"""
rule index_bam:
    input: "{Folder_name}/admixfrog_bams/{snpset}/{sample}.bam"
    output: "{Folder_name}/admixfrog_bams/{snpset}/{sample}.bam.bai"
    shell: "samtools index {input}"

def _prepare_bam_filter(wc):
    bam_path = config["EMH"][wc.sample]["bam"]
    s = "{}".format(bam_path)
    return s

def _unprocessed_bams_merge(wc):
    bam = config["EMH_unprocessed_bams"][wc.sample]["bam"]
    if "CHROM" in bam:
        bam = expand(bam, CHROM=[str(i+1) for i in range(22)] + ["X"])
    return bam

rule unprocessed_bams_merge:
    input:
        _unprocessed_bams_merge
    output:
        "{Folder_name}/admixfrog_bams/unprocessed_bams/{sample}.bam"
    run:
        if len(input) > 1:
            shell("samtools merge -c {output} {input}")
        else:
            shell("ln -sf {input} {output}")
"""
rule prepare_bam:
    input:
        bam= _prepare_bam_filter,
        bed= "/home/leonardo_iasi/EMH_Introgression_Project/Introgression_Detection/admixfrogREF/bed/{snpset}.bed",
    params:
        minlength=35,
        qual=25,
    output:
        bam = "{Folder_name}/admixfrog_bams/{snpset}/{sample}.bam",
        tbam=temp("{Folder_name}/bams_temp/{snpset}/{sample}.bam"),
    run:
        s = 'samtools view -h {input.bam} -L {input.bed} '
        s += "-q {params.qual} -F1 " #-b -o {output.bam}"
        s += "| awk ' {{if($1~/^@/ || length($10) >= {params.minlength}) print $0}}'"
        s += '|  samtools view - -S -h -b -o {output.tbam}'
        shell(s)
        shell("bam-rmdup -o {output.bam} {output.tbam}")

"""
# assumes bams are already deduplicated
rule prepare_bam_present_day:
    input:
        bam= _prepare_bam_filter,
        bed= "/home/leonardo_iasi/EMH_Introgression_Project/Introgression_Detection/admixfrogREF/bed/{snpset}.bed",
    params:
        minlength=35,
        qual=25,
    output:
        bam = "{Folder_name}/admixfrog_bams/{snpset}/{sample}.bam"
    run:
        s = 'samtools view -h {input.bam} -L {input.bed} '
        s += "-q {params.qual} "
        s += "| awk '{{ if ($0 ~ /^@/ || length($10) >= {params.minlength}) print $0}}'"
        s += '|  samtools view - -S -h -b -o {output.bam}'
        shell(s)


""" STATS """
rule compare_frag_call:
    input:
        true_rle =  "{Folder_name}/rle/{TAFP}/length_bin_{length_bin}/{bin_size}/{states}/{Tmap}/{sample}_{Tasc}.rle{TruePenalty}_downsampled.{Tdown}.xz",
        obs_rle =   "{Folder_name}/rle/{admf_params}/length_bin_{length_bin}/{bin_size}/{states}/{map}/{sample}_{snpset}.rle{penalty}_downsampled.{downsampled}.xz",
        #script_ = "fragment_calling_comparison.R"
    params:
        min_len = 0.0,
        type = "state",
        state = "NEA"
    output:
        res_obs = ("{Folder_name}/stats/Testing/{admf_params}/length_bin_{length_bin}/{bin_size}/{states}/TrueAFparams_{TAFP}/TrueAsc_{Tasc}/TrueDownsample_{Tdown}/True_Penalty_{TruePenalty}_minLen_0/TrueMap{Tmap}/{sample}_{snpset}.rle{penalty}_downsampled.{downsampled}_map{map}_frag_comp.csv")
    script:
        "fragment_calling_comparison.R"

def _inputfiles_plot_chromosome_comparison(wc):
    l = expand("{Folder_name}/stats/Testing/{admf_params}/length_bin_{length_bin}/{bin_size}/{states}/TrueAFparams_{TAFP}/TrueAsc_{Tasc}/TrueDownsample_{Tdown}/True_Penalty_{TruePenalty}_minLen_0/TrueMap{Tmap}/{sample}_{snpset}.rle{penalty}_downsampled.{downsampled}_map{map}_frag_comp.csv",
                Folder_name=Folder_name,admf_params=wc.admf_params,snpset=wc.snpset,sample = wc.sample,bin_size=wc.bin_size,length_bin=wc.length_bin,states=wc.states,TruePenalty=wc.TruePenalty,penalty=wc.penalty,downsampled=[1,0.2,0.1,0.02,0.01,0.005],Tasc=wc.Tasc,TAFP=wc.TAFP,Tdown=wc.Tdown,map=wc.map,Tmap=wc.Tmap)
    return l

rule plot_chromosome_comparison:
    input:
        true_rle =  expand("{Folder_name}/rle/{{TAFP}}/length_bin_{{length_bin}}/{{bin_size}}/{{states}}/{{Tmap}}/{{sample}}_{{Tasc}}.rle{{TruePenalty}}_downsampled.{{Tdown}}.xz",Folder_name=Folder_name),
        true_comp =   expand("{Folder_name}/stats/Testing/{{TAFP}}/length_bin_{{length_bin}}/{{bin_size}}/{{states}}/TrueAFparams_{{TAFP}}/TrueAsc_{{Tasc}}/TrueDownsample_{{Tdown}}/True_Penalty_{{TruePenalty}}_minLen_0/TrueMap{{Tmap}}/{{sample}}_{{Tasc}}.rle{{TruePenalty}}_downsampled.1_map{{map}}_frag_comp.csv",Folder_name=Folder_name),
        obs_comp =   _inputfiles_plot_chromosome_comparison,
        script_ = "plotting/plot_TechVal_chr_comparison.R"
    params:
        min_len = 0.0,
        min_len_obs = 0.0,
        type = "state",
        state = "NEA",
        chr = 2
    output:
        p_comp_table_data = temp("{Result_Folder_name}/CompFigures/{admf_params}/length_bin_{length_bin}/{bin_size}/{states}/TrueAFparams_{TAFP}/TrueAsc_{Tasc}/TrueDownsample_{Tdown}/True_Penalty_{TruePenalty}_minLen_0/TrueMap{Tmap}/{sample}_{snpset}.rle{penalty}_map{map}_FragCompTable.csv"),
    script:
        "plotting/plot_TechVal_chr_comparison.R"



def _inputfiles_plot_CLS_per_PenaltyParam(wc):
    l = expand("{Result_Folder_name}/CompFigures/{admf_params}/length_bin_{length_bin}/{bin_size}/{states}/TrueAFparams_{TAFP}/TrueAsc_{Tasc}/TrueDownsample_{Tdown}/True_Penalty_{TruePenalty}_minLen_0/TrueMap{Tmap}/{sample}_{snpset}.rle{penalty}_map{map}_FragCompTable.csv",
                Result_Folder_name=Result_Folder_name,admf_params=wc.admf_params,snpset=["ArchaicPlusUnfiltered",'DiagnosticSites','A1240k','archaicadmixtureAPX'],sample = wc.sample,bin_size=wc.bin_size,length_bin=wc.length_bin,states=wc.states,TruePenalty=wc.TruePenalty,penalty=[0.25],Tasc=wc.Tasc,TAFP=wc.TAFP,Tdown=wc.Tdown,map=wc.map,Tmap=wc.Tmap)
    return l

rule plot_CLS_per_PenaltyParam:
    input:
        p_comp_table_data = _inputfiles_plot_CLS_per_PenaltyParam,
        script_ = "plotting/plot_TechVal_freq_CLS_perPenaltyParam.R"
    output:
        p_frag_stats = "{Result_Folder_name}/CompFigures/{admf_params}/length_bin_{length_bin}/{bin_size}/{states}/TrueAFparams_{TAFP}/TrueAsc_{Tasc}/TrueDownsample_{Tdown}/True_Penalty_{TruePenalty}_minLen_0/TrueMap{Tmap}/PenaltyParam/{map}/{sample}_frag_comp_PenaltyParam.png",
        csv_frag_stats = "{Result_Folder_name}/CompFigures/{admf_params}/length_bin_{length_bin}/{bin_size}/{states}/TrueAFparams_{TAFP}/TrueAsc_{Tasc}/TrueDownsample_{Tdown}/True_Penalty_{TruePenalty}_minLen_0/TrueMap{Tmap}/PenaltyParam/{map}/{sample}_frag_comp_PenaltyParam.csv"
    script:
        "plotting/plot_TechVal_freq_CLS_perPenaltyParam.R"

def _inputfiles_stack_plot_CLS_per_frag_len(wc):
    l = expand("{Result_Folder_name}/CompFigures/{admf_params}/length_bin_{length_bin}/{bin_size}/{states}/TrueAFparams_{TAFP}/TrueAsc_{Tasc}/TrueDownsample_{Tdown}/True_Penalty_{TruePenalty}_minLen_0/TrueMap{Tmap}/PenaltyParam/{map}/{sample}_frag_comp_PenaltyParam.csv",
                Result_Folder_name=Result_Folder_name,admf_params=wc.admf_params,sample = wc.sample,bin_size=5000,length_bin=15,states="AFR-NEA-DEN",TruePenalty=wc.TruePenalty,Tasc=wc.Tasc,TAFP=wc.TAFP,Tdown=[1],map=["Shared_Map","AA_Map"],Tmap=["Shared_Map","AA_Map"])
    return l

rule stack_plot_CLS_per_frag_len:
    input:
        p_comp_table_data = _inputfiles_stack_plot_CLS_per_frag_len,
        meta_data = "/mnt/diversity/leonardo_iasi/EMH_Introgression_Project/New_REF_AA_with_APX_EMH_admixfrog0.7/stats/error2/AFR-NEA-DEN/QC_CurratedEMHPublished_length_bin_15_5000_archaicadmixtureAPX.wide_updated.csv",
        script_ = "plotting/stack_plot_CLS_per_frag_len.R"
    params:
        minLen = 0,
        maxLen = 2,
        plot_var = "Obs_Ascertainment",
        plot_asc = ""
    output:
        #p_frag_stats = "{Result_Folder_name}/CompFigures/CovAscParam{admf_params}/{map}/{sample}_TruePenalty_{TruePenalty}_TrueAsc_{Tasc}_TrueAfp_{TAFP}_TrueMap{Tmap}_ObsPenalty_{penalty}.png",
        p_frag_stats_RData = "{Result_Folder_name}/CompFigures/CovAscParam{admf_params}/{map}/Updated/{sample}_TruePenalty_{TruePenalty}_TrueAsc_{Tasc}_TrueAfp_{TAFP}_TrueMap{Tmap}_ObsPenalty_{penalty}.RData"
    script:
        "plotting/stack_plot_CLS_per_frag_len.R"


def exec_admixfrog(wc):
    Admixfrog_parameter = ["error2"]
    BinSize = 5000
    States = "-".join(["AFR", "NEA","DEN"])
    Sample = config['panels'][wc.panel]
    Snpset = ["ArchaicPlusUnfiltered",'DiagnosticSites','A1240k','archaicadmixtureAPX']
    penalty = [0.25]
    length_bin = 15
    downsampled = [1,0.2,0.1,0.02,0.01,0.005]
    Map = ["Shared_Map","AA_Map"]
    TMap = ["Shared_Map","AA_Map"]


    file_5 = expand("{Result_Folder_name}/CompFigures/CovAscParam{AF}/{M}/Updated/{C}_TruePenalty_{TP}_TrueAsc_{Tasc}_TrueAfp_{TAFP}_TrueMap{TM}_ObsPenalty_{E}.RData",
    Result_Folder_name=Result_Folder_name,C=Sample,AF=['error2'],TP=[0.25],E=penalty,Tasc=['DiagnosticSites','ArchaicPlusUnfiltered'],TAFP=["error2"],M=["Shared_Map","AA_Map"],TM=["Shared_Map","AA_Map"])

    


    files = file_5
    return files

rule run_series_cfg:
    input:
        exec_admixfrog
    output: touch('../controller/{panel}')
