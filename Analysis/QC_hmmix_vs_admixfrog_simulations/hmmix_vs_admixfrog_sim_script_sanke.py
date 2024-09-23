import os
import sys
import ts_helper_functions as ts_help
import demography_EMH_Introgression_Project_simple as EMH_demo
import msprime
import tskit
import pandas as pd
import numpy as np
from collections import defaultdict, namedtuple, Counter
from os import path
import gzip

configfile : "hmmix_vs_admixfrog_sim.yaml"

#################################################################################################
"""
Specify the number of diploid samples (N_target/2),
the number of replications and chromosomes here.
Additionaly, give the path to the folder th files can be stored.
"""
N_SAMPLES =  10
N_REPS = 10
CHROMS = [str(i+1) if i < 22  else "X" for i in range(0,22)]
Folder_name = '/mnt/diversity/leonardo_iasi/EMH_Introgression_Project/Simulations/RevisionSimulationswithDEN'
#################################################################################################



rule run_sim:
    input:
        #rec=get_rec_true
    priority : 1
    output:
        ts="{Folder_name}/sims/{sim}/{chrom}/{demo}.{rep}_ts.trees",
        target_sample_list="{Folder_name}/sims/{sim}/{chrom}/target_sample_list{demo}.{rep}.txt",
    run:
        demo = config['demography']['__default__'].copy()
        demo.update(config['series']['demography'][wildcards.demo])
        sim = config['sim']['__default__'].copy()
        sim.update(config['series']['sim'][wildcards.sim])

        map_ = sim['rec']

        if map_ is None:
            print("constant recombination used as true")
            rec =  None
        else:
            print("{0} empirical recombination map used as true".format(map_))
            rec = f"../Recombination_Maps/{map_}/chr{wildcards.chrom}.rec.gz"
        print(rec)
        #Altai + Chag + Vindija
        N_Na = 2
        s_Na = [125]
        N_Nc = 2
        s_Nc =  [90]
        N_Nv = 2
        s_Nv =  [55]
        sample_nea = N_Na * s_Na + N_Nc * s_Nc + N_Nv * s_Nv

        #Denisova 3
        N_D = 2
        s_D = [50]
        sample_den = N_D * s_D

        #Humans
        N_Hafr = demo['N_Hafr']
        N_Heur = demo['N_Heur']
        N_Hasn = demo['N_Hasn']
        sample_afr = [0] * N_Hafr
        sample_eur = [0] * N_Heur
        sample_asn = [0] * N_Hasn

        # target
        s_target =  demo['s_target']
        N_target = len(s_target) * 2
        p_target = demo['p_target']

        sample_target = np.repeat(s_target,N_target/len(s_target)).tolist()
        pop_target = np.repeat(p_target,N_target/len(p_target)).tolist()
        target_df = pd.DataFrame({'sample_age': sample_target, 'pop_target': pop_target})

        generation_time = demo['generation_time']
        admixture_proportion = demo['admixture_proportion']

        ind_pop_name = []
        sim_index = []
        if 'EUR' in target_df["pop_target"].values:
            print("Sampling ancient EUR")
            sample_eur = sample_eur + list(target_df[target_df.pop_target.isin(['EUR'])]['sample_age'])
            ind_pop_name.extend([f"EUR{i}" for i in range(N_Heur,(N_Heur + p_target.count('EUR') * 2))])
            sim_index.extend(range((len(sample_nea) + len(sample_den) + len(sample_afr) + N_Heur),
                              ((len(sample_nea) + len(sample_den) + len(sample_afr) + N_Heur) + p_target.count('EUR') * 2)))
        if 'ASN' in target_df["pop_target"].values:
            print("Sampling ancient ASN")
            sample_asn = sample_asn + list(target_df[target_df.pop_target.isin(['ASN'])]['sample_age'])
            ind_pop_name.extend([f"ASN{i}" for i in range(N_Hasn,(N_Hasn + p_target.count('ASN') * 2))])
            sim_index.extend(range((len(sample_nea) + len(sample_den) + len(sample_afr) + N_Heur + N_Hasn),
                              ((len(sample_nea) + len(sample_den) + len(sample_afr) + N_Heur + N_Hasn) + p_target.count('ASN') * 2)))

        target_df['ind_pop_name'] = ind_pop_name
        target_df['sim_index'] = sim_index
        target_df['haplo_name'] = [f"target{i}" for i in range(len(target_df['sim_index']))]
        target_df['sample_name'] = [f"target{i}" for i in np.repeat(range(N_SAMPLES),2)]
        target_df['true_admixture_time_in_gen'] = 2000
        target_df['sample_age_in_gen'] = target_df["sample_age"] * 1000 / generation_time
        target_df['generation_time'] = generation_time
        target_df['admixture_proportion'] = admixture_proportion


        #Contaminants
        N_cont = demo['N_cont']

        sample_list = {
              "NEA":  {"t_sample": sample_nea},
              "DEN":  {"t_sample": sample_den},
              "AFR": {"t_sample": sample_afr},
              "EUR": {"t_sample": sample_eur},
              "ASN": {"t_sample": sample_asn},
              "cont": {"t_sample": N_cont * [0]},
              "PAN": {"t_sample": 1 * [0]}
          }



        samples = EMH_demo.define_samples(sample_list,generation_time)

        #Demography= EMH_demo.demo_EMH_Introgression_Project(generation_time,admixture_proportion,True)
        #Demography= EMH_demo.demo_EMH_Introgression_Project_simplified(generation_time,admixture_proportion,True)
        admixture_proportion_den = 0.02
        Demography= EMH_demo.demo_EMH_Introgression_Project_with_diverged_DEN(generation_time,admixture_proportion,admixture_proportion_den,True)
        # the genome wide mutation rate
        mut_rate = sim['mut_rate']
        # the recombination rate can be either genome wide or a genetic map in HapMap format
        if rec is None:
            recomb_rate = sim['recomb_rate']
            target_df['recomb_rate'] = recomb_rate
            # the physical length of the chromosome in bp. If a genetic map is used remove this parameter
            # the length is then given by the length of the genetic map.
            seq_len = sim['seq_len']
            print("Simulation started")
            simulation_result_tree = msprime.sim_ancestry(
            samples=samples,
            sequence_length=seq_len,
            recombination_rate=recomb_rate,
            demography=Demography,
            record_migrations=True)
        else:
            print("Simulation started using empirical recombination map {0} for chr {1}".format(map_,wildcards.chrom))
            simulation_result_tree = msprime.sim_ancestry(
            samples=samples,
            recombination_rate=msprime.RateMap.read_hapmap(rec),
            demography=Demography,
            record_migrations=True)

        # Add mutations on tree sequence
        simulation_result = msprime.sim_mutations(simulation_result_tree, rate=mut_rate,model=msprime.BinaryMutationModel())
        print("Simulation finished")
        target_df.to_csv('{}'.format(output.target_sample_list), float_format="%.5f", index=False)
        simulation_result.dump('{}'.format(output[0]))


rule get_snp_table:
    input:
        ts="{Folder_name}/sims/{sim}/{chrom}/{demo}.{rep}_ts.trees",
        target_sample_list="{Folder_name}/sims/{sim}/{chrom}/target_sample_list{demo}.{rep}.txt"
    output:
        data="{Folder_name}/sims/{sim}/{chrom}/{demo}.{rep}_snps.tsv.gz"
    run:
        file_name = '{}'.format(output.data)
        ts=tskit.load('{}'.format(input.ts))
        target_sample_list = pd.read_csv('{}'.format(input.target_sample_list), sep=",")
        # Choose sample names (in correct order)
        sample_names = ["NEA","DEN","AFR","EUR","ASN","cont","PAN"]
        # Choose chromosome name
        chrom_name = wildcards.chrom
        target_df_d = dict(zip(target_sample_list.ind_pop_name,target_sample_list.haplo_name))
        ts_help.write_all_snps(
                        file_name=file_name,
                        ts=ts,
                        sample_names=sample_names,
                        chang_dict=target_df_d,
                        chrom=chrom_name)

rule get_true_fragments_simple:
    input:
        ts="{Folder_name}/sims/{sim}/{chrom}/{demo}.{rep}_ts.trees",
        target_sample_list="{Folder_name}/sims/{sim}/{chrom}/target_sample_list{demo}.{rep}.txt"
    output:
        TrueSeg=temp("{Folder_name}/sims/{sim}/{chrom}/TrueSegments_{demo}.{rep}.bed")
        #frags=temp("{Folder_name}/sims/{sim}/{chrom}/{demo}.{rep}_all_haplotypes.csv.gz")
    run:
        demo = config['demography']['__default__'].copy()
        demo.update(config['series']['demography'][wildcards.demo])
        target_sample_list = pd.read_csv('{}'.format(input.target_sample_list), sep=",")
        ts=tskit.load('{}'.format(input.ts))
        generation_time = demo['generation_time']
        NEA_segs = ts_help.get_introgressed_segments_new(ts,3,0,70000/int(generation_time),TestpopulationIDs=list(target_sample_list[target_sample_list.pop_target.isin(['EUR'])]['sim_index']))
        DEN_segs = ts_help.get_introgressed_segments_new(ts,3,7,49000/int(generation_time),TestpopulationIDs=list(target_sample_list[target_sample_list.pop_target.isin(['EUR'])]['sim_index']))
        with open('{}'.format(output.TrueSeg),'w') as out_true:
            for haplotype, archaic_segments in enumerate(NEA_segs):
                for archaic_segment in archaic_segments:
                    out_true.write('{}\tchr{}\t{}\t{}\tNEA\n'.format(haplotype,int(wildcards.chrom), int(archaic_segment[0]), int(archaic_segment[1])))
            for haplotype, archaic_segments in enumerate(DEN_segs):
                for archaic_segment in archaic_segments:
                    out_true.write('{}\tchr{}\t{}\t{}\tDEN\n'.format(haplotype,int(wildcards.chrom), int(archaic_segment[0]), int(archaic_segment[1])))
        

def merge_True_Archaic_Segments_files(wc):
    files = expand("{{Folder_name}}/sims/{{sim}}/{chrom}/TrueSegments_{{demo}}.{{rep}}.bed", chrom=CHROMS)
    return files

rule merge_True_Archaic_Segments_files:
    input:
        merge_True_Archaic_Segments_files
    output:
        temp("{Folder_name}/sims/{sim}/All_TrueSegments_{demo}.{rep}.bed")
    run:
        with open('{}'.format(output[0]),'w') as out_True_Archaic_Segments:
            for i, file in enumerate(input):
                with open(file) as data:
                    for line in data:
                        out_True_Archaic_Segments.write(line)

rule sample_info:
    input:
        "{Folder_name}/sims/{sim}/1/target_sample_list{demo}.{rep}.txt"
    output:
        "{Folder_name}/Sample_info{sim}.{demo}.{rep}.txt"
    run:
        s = "cat {input} > {output}"
        shell(s)

rule annotate_True_Archaic_Segments_files:
    input:
        segs="{Folder_name}/sims/{sim}/All_TrueSegments_{demo}.{rep}.bed",
        target_sample_list="{Folder_name}/sims/{sim}/1/target_sample_list{demo}.{rep}.txt"
    output:
        "{Folder_name}/Truth/{sim}/All_TrueSegmentsAnnotated_{demo}.{rep}.csv"
    run:
        target_sample_list = pd.read_csv('{}'.format(input.target_sample_list), sep=",")
        target_sample_list['haplo'] = target_sample_list.index
        segs = pd.read_csv('{}'.format(input.segs), sep="\t", names = ["haplo","chrom","start","end","source"])
        segs_anno = segs.merge(target_sample_list, on = "haplo", how = "left")
        segs_anno.to_csv('{}'.format(output), float_format="%.5f", index=False)

"""Admixfrog"""

rule admixfrog_input_ref:
    """create admixfrog input and panel file for 1 run 1 chromosome"""
    input:
        snps="{Folder_name}/sims/{sim}/{chrom}/{demo}.{rep}_snps.tsv.gz",
        #rec=get_rec_obs,
    resources:
        io=1
    priority : 0
    output:
        snps_asc=temp("{Folder_name}/sims/{sim}/{chrom}/{demo}.{seq_sim}.{rep}_snps_asc.tsv.gz"),
        ref=temp("{Folder_name}/sims/{sim}/{chrom}/{demo}.{seq_sim}.{rep}.panel.txt")

    run:
        seq_sim = config['seq_sim']['__default__'].copy()
        seq_sim.update(config['series']['seq_sim'][wildcards.seq_sim])
        sim = config['sim']['__default__'].copy()
        sim.update(config['series']['sim'][wildcards.sim])

        map_ = seq_sim['rec_obs']

        if map_ is None:
            print("constant recombination used as observed")
            rec_obs = None
        else:
            print("{0} chr {1} empirical recombination map used as observed".format(map_,wildcards.chrom))
            rec_obs = f"../Recombination_Maps/{map_}/chr{wildcards.chrom}.rec.gz"

        snps_table = pd.read_csv('{}'.format(input.snps), sep="\t")
        #scaled_map = pd.read_csv("Example_genetic_map_hapmap_formated.txt", sep="\t")

        # Choose sample names (in correct order)
        sample_names = ["NEA","DEN","AFR","EUR","ASN","target","cont","PAN"]

        chrom = wildcards.chrom

        if rec_obs is None:
            allele_table = ts_help.get_all_pop_allele_counts(
                                                  snps=snps_table,
                                                  sample_names=sample_names,
                                                  rec_rate=float(sim['recomb_rate'])*100,
                                                  rec_map=None,
                                                  chrom=chrom)
        else:
            rec_obs_map = pd.read_csv(rec_obs, sep="\t")
            allele_table = ts_help.get_all_pop_allele_counts(
                                                  snps=snps_table,
                                                  sample_names=sample_names,
                                                  rec_rate=None,
                                                  rec_map=rec_obs_map,
                                                  chrom=chrom)

        ascertainment = ts_help.ascertain_all_allele_counts(
                                                allele_table=allele_table,
                                                fixed_der = None if not seq_sim['fixed_der'] else seq_sim['fixed_der'],
                                                fixed_anc = None if not seq_sim['fixed_anc'] else seq_sim['fixed_anc'],
                                                polymorphic = None if not seq_sim['polymorphic'] else seq_sim['polymorphic'])

        ascertained_allele_table = ascertainment[0]
        ascertained_snp_table = snps_table[snps_table.index.isin(ascertainment[1])]
        ascertained_allele_table.columns = ascertained_allele_table.columns.str.replace(r"anc", "ref")
        ascertained_allele_table.columns = ascertained_allele_table.columns.str.replace(r"der", "alt")
        ascertained_allele_table.to_csv('{}'.format(output.ref), float_format="%.5f", index=False)
        ascertained_snp_table.to_csv('{}'.format(output.snps_asc), float_format="%.5f", index=False)



rule admixfrog_input_sample:
    """create admixfrog input and panel file for 1 run 1 chromosome"""
    input:
        snps_asc="{Folder_name}/sims/{sim}/{chrom}/{demo}.{seq_sim}.{rep}_snps_asc.tsv.gz",
        ref="{Folder_name}/sims/{sim}/{chrom}/{demo}.{seq_sim}.{rep}.panel.txt",
        target_sample_list="{Folder_name}/sims/{sim}/{chrom}/target_sample_list{demo}.{rep}.txt"
        #rec=get_rec_obs,
    resources:
        io=1
    priority : 0
    output:
        samples=temp("{Folder_name}/sims/{sim}/{chrom}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.sample.txt")


    run:
        demo = config['demography']['__default__'].copy()
        demo.update(config['series']['demography'][wildcards.demo])

        ascertained_snp_table = pd.read_csv('{}'.format(input.snps_asc), sep=",")
        ascertained_allele_table = pd.read_csv('{}'.format(input.ref), sep=",")

        # Choose sample names (in correct order)
        sample_names = ["NEA","DEN","AFR","EUR","ASN","target","cont","PAN"]

        target_sample_list = pd.read_csv('{}'.format(input.target_sample_list), sep=",")

        chrom = wildcards.chrom

        ID = wildcards.id

        ascertained_allele_table.columns = ascertained_allele_table.columns.str.replace(r"ref", "anc")
        ascertained_allele_table.columns = ascertained_allele_table.columns.str.replace(r"alt", "der")

        samples = range(demo['N_target'])
        n_samples = len(samples) // 2
        samples = np.array(samples).reshape(n_samples, 2).tolist()
        ids = [f'target{i}' for i in samples[int(ID)]]
        print("Generating inputfiles for gtmode")
        sample_input=ts_help.admixfrog_input_gt(
                                ids=ids,
                                allele_table=ascertained_allele_table,
                                snp_table=ascertained_snp_table)
        sample_input.to_csv('{}'.format(output.samples), float_format="%.5f", index=False)


rule merge_ref:
    input:
        expand("{{Folder_name}}/sims/{{sim}}/{chrom}/{{demo}}.{{seq_sim}}.{{rep}}.panel.txt", chrom=CHROMS)
    priority : 100
    group : "merge"
    output:
        "{Folder_name}/infiles/{sim}/{demo}.{seq_sim}.{rep}.ref.xz"
    shell:
        "head -qn1 {input[0]} | xz -c > {output} && "
        "tail -qn+2 {input} | xz -c >> {output} "

rule merge_sample:
    input:
        #expand("{{Folder_name}}/sims/{{sim}}/{chrom}/{{demo}}.{{seq_sim}}.{{rep}}.{{id}}.sample.txt", chrom=CHROMS)
        expand("{{Folder_name}}/sims/{{sim}}/{chrom}/admixforg_params{{admixP}}/{{demo}}.{{seq_sim}}.{{rep}}.{{id}}.sample.txt", chrom=CHROMS)
    priority : 1020
    group : "merge"
    output:
        "{Folder_name}/infiles/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.in.xz"
    shell:
        "head -qn1 {input[0]} | xz -c > {output} && "
        "tail -qn+2 {input} | xz -c >> {output} "

def get_recs(wc):
    pars = config['sim']['__default__'].copy()
    pars.update(config['series']['sim'][wc.sim])
    return expand("/mnt/diversity/bpeter/admixfrog/sims3/recs/hapmap/{map}/chr{chrom}.rec.gz",
                 map=pars['rec'],
                 chrom = CHROMS)

rule merge_true_reps:
    input:
        f=expand("{{Folder_name}}/sims/{{sim}}/{{demo}}.{rep}_all_haplotypes.xz",
            rep=range(N_REPS)),
        idtbl ='{Folder_name}/sims/{sim}/{demo}.idtbl'
    output:
        true_frags="{Folder_name}/sims/true/{sim}/{demo}_all_haplotypes_merged.xz"
    script: 'scripts/merge_true2.R'

rule run_admixfrog_gtmode:
    input:
        infile="{Folder_name}/infiles/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.in.xz",
        ref="{Folder_name}/infiles/{sim}/{demo}.{seq_sim}.{rep}.ref.xz"
    params:
        ll_tol = 1e-2,
        freq_f = 3,
        freq_c = 3,
        error = 1e-2,
        ancestral = "PAN",
        max_iter = 250,
        n_post_rep = 200
    log :
        log="{Folder_name}/admixfrog/{sim}/gtmode/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.log"
    output:
        res_bin = "{Folder_name}/admixfrog/gtmode/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.bin.xz",
        res_snp = "{Folder_name}/admixfrog/gtmode/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.snp.xz",
        res_cont= "{Folder_name}/admixfrog/gtmode/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.cont.xz",
        res_par = "{Folder_name}/admixfrog/gtmode/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.pars.yaml",
        res_rle = "{Folder_name}/admixfrog/gtmode/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.rle.xz",
        res_res = "{Folder_name}/admixfrog/gtmode/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.res.xz",
        res2_res = "{Folder_name}/admixfrog/gtmode/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.res2.xz",
    run:
        admixfrog = config['admixfrog']['__default__'].copy()
        admixfrog.update(config['series']['admixfrog'][wildcards.admixP])
        outname = path.splitext(path.splitext(output.res_bin)[0])[0]
        #states = wildcards.states.split("-")
        states = admixfrog['state_ids']
        cont_id = admixfrog['cont_id']
        run_penalty = admixfrog['run_penalty']
        bin_size = 5000
        n_CHR = len(CHROMS)

        s = "admixfrog --infile {input.infile} --ref {input.ref} -o {outname} "
        s += "--male --filter-pos 50 --filter-map 0"
        s += " --gt-mode "
        s += " --states {states} "
        s += " --ll-tol {params.ll_tol} "
        s += " --bin-size {bin_size} "
        s += " --est-F --est-tau --freq-F {params.freq_f} "
        s += " --dont-est-contamination "
        s += " --c0 0"
        s += " --e0 {params.error} "
        s += " --ancestral {params.ancestral} "
        s += " --run-penalty {run_penalty}"
        s += " --max-iter {params.max_iter}"
        s += " --n-post-replicates {params.n_post_rep}"
        s += " --init-guess NEA"
        s += " -P" # position mode
        if n_CHR < 23:
            s += " --chroms {n_CHR}"
        shell(s)

rule call_runs:
    input:
        rle = "{Folder_name}/admixfrog/gtmode/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.bin.xz"

    output:
        rle = "{Folder_name}/rle/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.rle{rle}.xz",

    run:
        run_penalty = wildcards.rle

        s = "admixfrog-rle --out {output.rle} --in {input.rle} "
        s += "--run-penalty {run_penalty} "
        shell(s)

"""hmmix"""

rule remove_outgroup:
    input:
        ts="{Folder_name}/sims/{sim}/{chrom}/{demo}.{rep}_ts.trees"
    output:
        All_Ind_obs="{Folder_name}/sims/{sim}/{chrom}/{demo}.{rep}_ts_Outgroup_removed.trees",
    run:
        ts=tskit.load('{}'.format(input.ts))
        Outgroup_removed = ts_help.remove_outgroup(ts,2)
        Outgroup_removed.dump('{}'.format(output[0]))

rule get_obs_file:
    input:
        ts="{Folder_name}/sims/{sim}/{chrom}/{demo}.{rep}_ts_Outgroup_removed.trees",
        target_sample_list="{Folder_name}/sims/{sim}/{chrom}/target_sample_list{demo}.{rep}.txt"
    output:
        All_Ind_obs=temp("{Folder_name}/sims/{sim}/{chrom}/Observations_{demo}.{rep}.txt.gz"),
    run:
        demo = config['demography']['__default__'].copy()
        demo.update(config['series']['demography'][wildcards.demo])
        Outgroup_removed = ts=tskit.load('{}'.format(input.ts))
        target_sample_list = pd.read_csv('{}'.format(input.target_sample_list), sep=",")    
        Testpopulation = target_sample_list['sim_index'].values   
        """ Caution Archaic sample ids hardcoded!!! """
        ALT = [0,1]
        CHA = [2,3]
        VIN = [4,5]
        DEN = [6,7]
        with gzip.open('{}'.format(output.All_Ind_obs),'wt') as out_obs:
            for variant in Outgroup_removed.variants():
                for Ind in range(len(Testpopulation)):
                    haplotypes_as_string = variant.genotypes[Testpopulation[int(Ind)]]
                    if haplotypes_as_string == 0:
                        pass
                    else:
                        archaicsnp = []
                        if variant.genotypes[Testpopulation[int(Ind)]] == any(variant.genotypes[ALT]):
                            archaicsnp += ["ALT"]
                        if variant.genotypes[Testpopulation[int(Ind)]] == any(variant.genotypes[CHA]):
                            archaicsnp += ["CHA"]
                        if variant.genotypes[Testpopulation[int(Ind)]] == any(variant.genotypes[VIN]):
                            archaicsnp += ["VIN"]
                        if variant.genotypes[Testpopulation[int(Ind)]] == any(variant.genotypes[DEN]):
                            archaicsnp += ["DEN"]
                        out_obs.write('{}\tchr{}\t{}\t{}\tA\tC\t{}\n'.format(int(Ind),wildcards.chrom,int(variant.site.position),int(variant.site.mutations[0].time),",".join(archaicsnp)))
                        #out_obs.write('{}\tchr{}\t{}\tA\tC\n'.format(int(Ind),wildcards.chrom,int(variant.site.position)))
                        #out_obs.write('Ind_{}\tchr{}\t{}\t{}\tA\tC\n'.format(int(Ind),wildcards.chrom,int(variant.site.position),int(variant.mutations.time)))

rule merge_all_sample_all_chrom:
    input:
        expand("{{Folder_name}}/sims/{{sim}}/{chrom}/Observations_{{demo}}.{{rep}}.txt.gz", chrom=CHROMS)
    output:
        "{Folder_name}/sims/{sim}/All_Observations_{demo}.{rep}.txt.gz"
    shell:
        "zcat {input} | gzip > {output}"

rule hmmix_obs_ind_input:
    input:
        "{Folder_name}/sims/{sim}/All_Observations_{demo}.{rep}.txt.gz"
    output:
        #temp("{Folder_name}/sims/{sim}/{demo}/{rep}/Obs/obs.Ind_{Ind}.txt"),
        Obs = "{Folder_name}/sims/{sim}/{demo}/{rep}/Obs/obs.Ind_{id}.txt",
        Anno_Obs = "{Folder_name}/sims/{sim}/{demo}/{rep}/Anno_Obs/obs.Ind_{id}.txt"
    run:
        Ind = [int(wildcards.id) * 2, int(wildcards.id) * 2 + 1]
        obs_file = pd.read_csv('{}'.format(input),sep="\t",names= ["Ind","chrom" , "pos" , "time", "ancestral_base" , "genotype", "archaicsnps"])
        obs_file_hap_1 = obs_file[(obs_file.Ind == Ind[0])]
        obs_file_hap_2 = obs_file[(obs_file.Ind == Ind[1])]

        df = pd.merge(obs_file_hap_1, obs_file_hap_2, on=['chrom','pos'], how='outer')
        df["ancestral_base"] = "A"
        df['genotype_x'].fillna(value="A", inplace=True)
        df['genotype_y'].fillna(value="A", inplace=True)
        df['time_x'].fillna(value=df['time_y'], inplace=True)
        df['archaicsnps_x'].fillna(value=df['archaicsnps_y'], inplace=True)
        df['genotype'] = [''.join(i) for i in zip(df['genotype_x'], df['genotype_y'])]
        df['time'] = df['time_x']
        df['archaicsnps'] = df['archaicsnps_x']
        res_df_anno = df[["chrom" , "pos" , "time", "ancestral_base" , "genotype", "archaicsnps"]]
        res_df_anno.to_csv("{}".format(output.Anno_Obs),  sep = "\t",header = False, index = False)
        res_df = res_df_anno[['chrom','pos','ancestral_base','genotype']]
        res_df.to_csv("{}".format(output.Obs),  sep = "\t",header = False, index = False)

rule train_hmmix:
    input:
        "{Folder_name}/sims/{sim}/{demo}/{rep}/Obs/obs.Ind_{id}.txt",
    output:
        temp("{Folder_name}/sims/{sim}/Traind_hmmix_Ind_{id}.{demo}.{rep}.json")

    run:
        s = "source activate hmmix; hmmix train -obs={0} -param=../Initialguesses.json -haploid -out={1}; conda deactivate".format(input[0],output[0])
        print(s)
        shell(s)

rule decode_hmmix:
    input:
        obs = "{Folder_name}/sims/{sim}/{demo}/{rep}/Obs/obs.Ind_{id}.txt",
        param = "{Folder_name}/sims/{sim}/Traind_hmmix_Ind_{id}.{demo}.{rep}.json"
    output:
        out_hap1 = "{Folder_name}/hmmix/{sim}/{demo}/{rep}/Decode/Ind_{id}.hap1.txt",
        out_hap2 = "{Folder_name}/hmmix/{sim}/{demo}/{rep}/Decode/Ind_{id}.hap2.txt"
    run:
        s = "source activate hmmix; hmmix decode -obs={0} -param={1} -haploid -out={2}; conda deactivate".format(input.obs,input.param,output[0][:-9])
        print(s)
        shell(s)

rule match_hmmix_calls:
    input:
        hmmix_hap1 =  "{Folder_name}/hmmix/{sim}/{demo}/{rep}/Decode/Ind_{id}.hap1.txt",
        hmmix_hap2 = "{Folder_name}/hmmix/{sim}/{demo}/{rep}/Decode/Ind_{id}.hap2.txt",
        Anno_Obs = "{Folder_name}/sims/{sim}/{demo}/{rep}/Anno_Obs/obs.Ind_{id}.txt"
    output:
        hmmix_hap1matched =  "{Folder_name}/hmmix/{sim}/{demo}/{rep}/Decode/Ind_{id}.cutoff{hmmix_cutoff}.matched.csv",
    script:
        "match_hmmix.R"

"""comparison"""

rule compare_frag_call_between:
    input:
        rle_admixfrog = "{Folder_name}/rle/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.rle{rle}.xz",
        hmmix_hap1matched =  "{Folder_name}/hmmix/{sim}/{demo}/{rep}/Decode/Ind_{id}.cutoff{hmmix_cutoff}.matched.csv",
    output:
        res_obs = "{Folder_name}/comparison/{sim}/{demo}/admixforg_params{admixP}/{seq_sim}/{source}/BetweenComp.rle{rle}.minlen{minlen}.cutoff{hmmix_cutoff}.matched{hmmix_matched}.{rep}.{id}.csv"
    script:
        "compare_hmmix_with_admixfrog.R"

def merge_all_res_files_between(wc):
    files = expand("{{Folder_name}}/comparison/{{sim}}/{{demo}}/admixforg_params{{admixP}}/{{seq_sim}}/{{source}}/BetweenComp.rle{{rle}}.minlen{{minlen}}.cutoff{{hmmix_cutoff}}.matched{{hmmix_matched}}.{rep}.{id}.csv",
                   rep=range(N_REPS),id=range(N_SAMPLES))
    return files

rule merge_frag_call_between:
    input:
        merge_all_res_files_between
    output:
        "{Folder_name}/comparison/{sim}/{demo}/admixforg_params{admixP}/{seq_sim}/{source}/BetweenComparison{name}.rle{rle}.minlen{minlen}.cutoff{hmmix_cutoff}.matched{hmmix_matched}frag_comp_res.csv"
    shell:
        "head -qn1 {input[0]}  > {output} && "
        "tail -qn+2 {input}  >> {output} "

rule compare_frag_call:
    input:
        true_rle =  "{Folder_name}/Truth/{sim}/All_TrueSegmentsAnnotated_{demo}.{rep}.csv",
        rle_admixfrog = "{Folder_name}/rle/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.rle{rle}.xz",
        hmmix_hap1matched =  "{Folder_name}/hmmix/{sim}/{demo}/{rep}/Decode/Ind_{id}.cutoff{hmmix_cutoff}.matched.csv",
    output:
        res_obs = "{Folder_name}/comparison/{sim}/{demo}/admixforg_params{admixP}/{seq_sim}/{source}/SegmentComparison.rle{rle}.minlen{minlen}.cutoff{hmmix_cutoff}.matched{hmmix_matched}.{rep}.{id}.csv"
    script:
        "compare_hmmix_admixfrog_segments.R"

def merge_all_res_files(wc):
    files = expand("{{Folder_name}}/comparison/{{sim}}/{{demo}}/admixforg_params{{admixP}}/{{seq_sim}}/{{source}}/SegmentComparison.rle{{rle}}.minlen{{minlen}}.cutoff{{hmmix_cutoff}}.matched{{hmmix_matched}}.{rep}.{id}.csv",
                   rep=range(N_REPS),id=range(N_SAMPLES))
    return files

rule merge_frag_call:
    input:
        merge_all_res_files
    output:
        "{Folder_name}/comparison/{sim}/{demo}/admixforg_params{admixP}/{seq_sim}/{source}/Comparison{name}.rle{rle}.minlen{minlen}.cutoff{hmmix_cutoff}.matched{hmmix_matched}frag_comp_res.csv"
    shell:
        "head -qn1 {input[0]}  > {output} && "
        "tail -qn+2 {input}  >> {output} "

def exec_sim(wc):
    sim = config['set'][wc.name]['sim']
    demo = config['set'][wc.name]['demography']
    seq_sim = config['set'][wc.name]['seq_sim']
    admixP = config['set'][wc.name]['admixfrog']
    rle = 0.25
    minlen = [0]
    hmmix_cutoff = [0.8]
    hmmix_matched = ["Matched","All"]
    source = ["NEA","DEN"]



    file_1 = expand("{Folder_name}/Truth/{sim}/All_TrueSegmentsAnnotated_{demo}.{rep}.csv",Folder_name=Folder_name,
    sim=sim,demo=demo,seq_sim=seq_sim,rep=range(N_REPS),id=range(N_SAMPLES))
    file_2 = expand("{Folder_name}/rle/{sim}/admixforg_params{admixP}/{demo}.{seq_sim}.{rep}.{id}.rle{rle}.xz",
    Folder_name=Folder_name,sim=sim,admixP=admixP,demo=demo,seq_sim=seq_sim,rep=range(N_REPS),id=range(N_SAMPLES),rle=rle)
    file_3 = expand("{Folder_name}/Sample_info{sim}.{demo}.{rep}.txt",Folder_name=Folder_name,
    sim=sim,demo=demo,rep=range(N_REPS))

    file_4 = expand("{Folder_name}/hmmix/{sim}/{demo}/{rep}/Decode/Ind_{id}.hap1.txt",
                    Folder_name=Folder_name,sim=sim,demo=demo,rep=range(N_REPS),id=range(N_SAMPLES))
    file_5 = expand("{Folder_name}/hmmix/{sim}/{demo}/{rep}/Decode/Ind_{id}.cutoff{hmmix_cutoff}.matched.csv",
                     Folder_name=Folder_name,sim=sim,demo=demo,rep=range(N_REPS),id=range(N_SAMPLES),hmmix_cutoff=hmmix_cutoff,)

    comp_file = expand("{Folder_name}/comparison/{sim}/{demo}/admixforg_params{admixP}/{seq_sim}/{source}/Comparison{name}.rle{rle}.minlen{minlen}.cutoff{hmmix_cutoff}.matched{hmmix_matched}frag_comp_res.csv",
                       Folder_name=Folder_name,sim=sim,admixP=admixP,demo=demo,seq_sim=seq_sim,rep=range(N_REPS),id=range(N_SAMPLES),source=source,rle=rle,minlen=minlen,hmmix_cutoff=hmmix_cutoff,hmmix_matched=hmmix_matched,name=wc.name)
    comp_file2 = expand("{Folder_name}/comparison/{sim}/{demo}/admixforg_params{admixP}/{seq_sim}/{source}/BetweenComparison{name}.rle{rle}.minlen{minlen}.cutoff{hmmix_cutoff}.matched{hmmix_matched}frag_comp_res.csv",
                       Folder_name=Folder_name,sim=sim,admixP=admixP,demo=demo,seq_sim=seq_sim,rep=range(N_REPS),id=range(N_SAMPLES),source=source,rle=rle,minlen=minlen,hmmix_cutoff=hmmix_cutoff,hmmix_matched=hmmix_matched,name=wc.name)

    files =  file_1 + file_2 + file_3 + file_4 + file_5 + comp_file

    files = file_1 + file_2 + file_5 + comp_file + comp_file2

    return files

rule run_series_cfg:
    input:
        exec_sim
    output:
        touch("../simulations/{name}")
