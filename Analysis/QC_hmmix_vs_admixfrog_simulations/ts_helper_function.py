"""
Copyright (c) 2021, Leonardo Iasi, Laurits Skov and Benjamin Peter (Max Planck Institute for evolutionary anthropology)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
All advertising materials mentioning features or use of this software must display the following acknowledgement: This product includes software developed by the copyright holders: Leonardo Iasi, Laurits Skov and Benjamin Peter (Max Planck Institute for evolutionary anthropology).
Neither the name of Leonardo Iasi, Laurits Skov and Benjamin Peter (Max Planck Institute for evolutionary anthropology) nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY Leonardo Iasi, Laurits Skov and Benjamin Peter (Max Planck Institute for evolutionary anthropology) AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL Leonardo Iasi, Laurits Skov and Benjamin Peter (Max Planck Institute for evolutionary anthropology) BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


import msprime
import gzip
import pandas as pd
from scipy.stats import binom, poisson
import numpy as np
from scipy import stats
import tskit
import numpy as np
import collections


def write_all_snps(file_name, ts, sample_names, chang_dict=None, chrom='1'):
    """
    This function writes genotype matrix to a gzipped file with user defined pop names.
    The population names must be given in order of msprime sample id (ascending order)
    """
    """Generates the individual sample names from the user definded pop names"""
    samples = []
    for p in enumerate(sample_names):
        n_pop = len(ts.get_samples(p[0]))
        samples.extend([f"{p[1]}{i}" for i in range(n_pop)])
    if chang_dict is not None:
        samples = [chang_dict.get(item,item)  for item in samples]
    """Write genotype matrix to file."""
    with gzip.open(file_name, 'wt') as f:
        f.write("chrom\tsite_id\tpos\t")
        f.write("\t".join(samples))
        f.write("\n")
        for v in ts.variants():
            gt = '\t'.join(str(gt) for gt in v.genotypes)
            f.write(f"{chrom}\t{v.site.id}\t{int(v.site.position)}\t{gt}\n")



def get_all_pop_allele_counts(snps, sample_names, rec_rate=None, rec_map=None, chrom='1'):
    """
    This function creates a panda table with allele counts per SNP for each population from a
    snps file table (created by write_all_snps function)
    The population names must be given in order of msprime sample id (ascending order)
    """

    """"Reading in data from csv file containing the variants (snp file)"""
    if 'pos' not in snps:
        snps['pos'] = snps.index.astype(int)
    snps.pos = snps.pos.astype(int)
    snps2 = snps.drop_duplicates(subset=['pos'])


    D = dict()
    D['chrom'] = chrom
    D['site_id'] = snps2.site_id
    D['pos'] = snps2.pos.astype(int)
    if all(i is None for i in [rec_map,rec_rate]):
        D['map'] = snps2.pos * 1e-8
    elif rec_map is None:
        D['map'] = snps2.pos * rec_rate
    else:
        D['map'] = np.interp(snps2.pos, rec_map.pos, rec_map.map)

    D['anc'] = 'A'
    D['der'] = 'G'

    for p in sample_names:
        pop_x_cols = [col for col in snps2.columns if col.startswith(f"{p}")]
        n_pop_x = len(pop_x_cols)
        """
        Determining the count of derived (1) and ancestral (0) alleles per population:
        the number of haploid individuals minus the number of derived alleles determines
        the number of ancestral alleles per population
        """
        D[f"{p}_der"] = np.sum(snps2[pop_x_cols], 1)
        D[f"{p}_anc"] = n_pop_x - D[f"{p}_der"]


    """
    Putting it into a data frame and check if number of input and poutput snps is the same
    and dropping potential duplicates
    """
    allele_table = pd.DataFrame.from_dict(D)
    assert allele_table.shape[0] == snps2.shape[0]
    assert allele_table.drop_duplicates(subset=['chrom', 'pos']).shape[0] == allele_table.shape[0]
    return(allele_table)


def ascertain_all_allele_counts(allele_table, fixed_der=None,fixed_anc=None,polymorphic=None,polarized=False):
    """
    This function ascertaines the allele table object (generated by the get_all_pop_allele_counts function)
    by user defined populations and outputs the filtered table.
    The populations ascertainment can be done using different categories:
    either be fixed for the derived/ancestral allele or polymorphic or both.
    Multiple populations can be defined per category.
    Optionally a snp file can be provided and is filtered in the same way and outputted.
    """

    """First the variants are filtered out which are not fixed"""
    if all(v is None for v in (fixed_der,fixed_anc)):
        filter_1 = (allele_table.pos == allele_table.pos)
    else:
        if fixed_der == None:
            fixed_anc_der = [f"{i}_der" for i in fixed_anc]
            fixed_anc_anc = [f"{i}_anc" for i in fixed_anc]
            print([f"Ascertaining for fixed ancestral sites in {i}" for i in fixed_anc])
            filter_1 = (allele_table[fixed_anc_der].sum(axis=1) == 0)
        if fixed_anc == None:
            fixed_der_der = [f"{i}_der" for i in fixed_der]
            fixed_der_anc = [f"{i}_anc" for i in fixed_der]
            print([f"Ascertaining for fixed derived sites in {i}" for i in fixed_der])
            filter_1 = (allele_table[fixed_der_anc].sum(axis=1) == 0)
        elif all(v is not None for v in (fixed_der,fixed_anc)):
            fixed_anc_der = [f"{i}_der" for i in fixed_anc]
            fixed_anc_anc = [f"{i}_anc" for i in fixed_anc]
            fixed_der_der = [f"{i}_der" for i in fixed_der]
            fixed_der_anc = [f"{i}_anc" for i in fixed_der]
            if polarized == True:
                print([f"Ascertaining for fixed derived sites in {i}" for i in fixed_der])
                print([f"Ascertaining for fixed ancestral sites in {i}" for i in fixed_anc])
                print("non polarized ascertainment")
                filter_1 = (allele_table[fixed_anc_der].sum(axis=1) + allele_table[fixed_der_anc].sum(axis=1) == 0)
            else:
                print([f"Ascertaining for fixed derived sites in {i}" for i in fixed_der])
                print([f"Ascertaining for fixed ancestral sites in {i}" for i in fixed_anc])
                print("polarized ascertainment")
                filter_1 = (allele_table[fixed_anc_der].sum(axis=1) + allele_table[fixed_der_anc].sum(axis=1) == 0) | (allele_table[fixed_der_der].sum(axis=1) + allele_table[fixed_anc_anc].sum(axis=1) == 0)
    allele_table_asc_1 = allele_table[filter_1]

    """Second, filtering for polymorphic variants"""
    if polymorphic == None:
        filter_2 = (allele_table_asc_1.pos == allele_table_asc_1.pos)
    else:
        polymorphic_der = [f"{i}_der" for i in polymorphic]
        polymorphic_anc = [f"{i}_anc" for i in polymorphic]
        filter_2 = (allele_table_asc_1[polymorphic_der].sum(axis=1) !=0) & (allele_table_asc_1[polymorphic_anc].sum(axis=1) !=0)
        print([f"Ascertaining for polymorphic sites in {i}" for i in polymorphic])
    allele_table_asc_2 = allele_table_asc_1[filter_2]
    filter_all=allele_table_asc_2.index
    if allele_table_asc_2.shape[0] == allele_table.shape[0]:
        print("Data not ascertained")
    return(allele_table_asc_2,filter_all)

def admixfrog_sample_in(sample_reads,min_len,bin_len):
    max_len= sample_reads.len
    sample_reads.loc[sample_reads['deam'] == 0, 'deam'] = "deam"
    sample_reads.loc[sample_reads['deam'] == -1, 'deam'] = "nodeam"
    i = 1
    while (min_len + bin_len*i) <= max(max_len + bin_len):
        sample_reads.loc[(sample_reads['len'] >= (min_len + bin_len*(i-1))) & (sample_reads['len'] <= (min_len + bin_len*(i))), 'len'] = i-1
        i=i+1
    admixfrog_sample_in = sample_reads[['chrom', 'pos']].copy()
    admixfrog_sample_in['lib'] = [p1 + '_' + p2 + '_' + p3 for p1, p2, p3 in zip(sample_reads['lib'], sample_reads['len'].astype(str), sample_reads['deam'])]
    admixfrog_sample_in['tref'] = sample_reads['tref']
    admixfrog_sample_in['talt'] = sample_reads['talt']
    admixfrog_sample_in['tdeam'] = np.repeat(0, admixfrog_sample_in.shape[0], axis=0)
    admixfrog_sample_in['tother'] = np.repeat(0, admixfrog_sample_in.shape[0], axis=0)
    admixfrog_sample_in=admixfrog_sample_in.groupby(['chrom','pos','lib']).sum().reset_index()
    return admixfrog_sample_in

def admixfrog_input_gt(ids, allele_table ,snp_table):
    """
    This function takes 1 or more diploid individuals specified by their haploid ids from the snp file,
    the allele table file and the snp file to simulates coverage for each SNP and the composition of
    endogenous and contaminant read per SNP observed. The user defines the source of contamination by
    its pop lable. The coverage, contamination_percent and libs (string) parameters must be given as arrays.
    """
    if isinstance(len(ids)/2, float) is not True:
        raise ValueError("Number of ids must be two or multiples of two")
    S = []

    data = allele_table[['chrom', 'pos']].copy()
    data['talt'] = np.sum(snp_table[ids],1)
    data['tref'] = 2 - data['talt']

    """Excludes all SNPs with 0 reads"""
    data = data[data.tref+data.talt>0]
    S.append(data)
    """Write the coverage file to csv"""
    data = pd.concat(S).sort_values(['chrom', 'pos'])

    return(data)

def get_introgressed_segments_new(ts,mixed_pop_ID,introgressing_pop_ID,split_time,TestpopulationIDs=None,mig_events_list=None,file_name=None,chrom='1'):
    """
    This function takes a tree sequence object (while running the msprim sim, record_migrations must have been set to True!!!),
    and returns all (merged) segments introgressed from the introgressing_pop found in sampled individuals of the mixed_pop.
    If a file_name is given it writes a file with individual id, chrom, start, end of intriogressed segment.
    """
    if TestpopulationIDs is not None:
        Testpopulation = TestpopulationIDs
        print("Taking subset individuals from Testpopulation")
    else:
        Testpopulation = ts.get_samples(mixed_pop_ID)
        print("Taking all sampled individuals from Testpopulation")

    de_seg = {i: [] for i in Testpopulation}

    """
    This is a helper function and merges neighboring segments in an individuals genome to one.
    """
    def combine_segs(segs, get_segs = False):
        merged = np.empty([0, 2])
        if len(segs) == 0:
            if get_segs:
                return([])
            else:
                return(0)
        sorted_segs = segs[np.argsort(segs[:, 0]), :]
        for higher in sorted_segs:
            if len(merged) == 0:
                merged = np.vstack([merged, higher])
            else:
                lower = merged[-1, :]
                if higher[0] <= lower[1]:
                    upper_bound = max(lower[1], higher[1])
                    merged[-1, :] = (lower[0], upper_bound)
                else:
                    merged = np.vstack([merged, higher])
        if get_segs:
            return(merged)
        else:
            return(np.sum(merged[:, 1] - merged[:, 0])/seq_len)
    """
    Looping through the tree sequence tables of migrants. If a migrant is from the
    introgressing_pop found in the mixed_pop, the sub interval of the tree is taken.
    From this interval all leaves from gthe migrant node (i.e. descendedts of the migrant)
    and are appended to our introgressed seg object if they are in any sampled ind of the mixed_pop.
    """
    mig_events = []
    for mr in ts.migrations():
        if mr.source == mixed_pop_ID and mr.dest == introgressing_pop_ID:
            if mr.time == split_time: # migraton caused by merge
                pass
            mig_events.append(mr)

    if mig_events_list is not None:
        print("Conditioning on passed migration events list")
        mig_events_c = []
        for mr in mig_events:
             if mr.node in mig_events_list:
                mig_events_c.append(mr)
        mig_events = mig_events_c

    mig_events.sort(key = lambda x: x.left)
    mig_events.sort(key = lambda x: x.right)
    for mr in mig_events:
        for tree in ts.keep_intervals([[mr.left,mr.right]],simplify=False,record_provenance=True).trees(leaf_lists=True):
            if mr.left > tree.get_interval()[0]:
                continue
            if mr.right <= tree.get_interval()[0]:
                break
            for l in tree.leaves(mr.node):
                if l in Testpopulation:
                    de_seg[l].append(tree.get_interval())

    """
    Merging all length neighboring intervals of introgressed segments using the combine_segs
    helper function.
    """
    true_de_segs = [combine_segs(np.array(de_seg[i]), True) for i in sorted(de_seg.keys())]
    """
    Writing the introgressed segs of a individual (haploid genome id) with chrom, and start and stop in bp to a file.
    """
    if file_name is not None:
        with open('{}'.format(file_name),'w') as out_true:
            for haplotype, archaic_segments in enumerate(true_de_segs):
                for archaic_segment in archaic_segments:
                    out_true.write('{}\tchr{}\t{}\t{}\n'.format(haplotype,int(chrom), int(archaic_segment[0]), int(archaic_segment[1])))
    else:
        return(true_de_segs)


def remove_outgroup(ts,outgroup_ID):
    """
    This function takes a tree sequence object and removes all variants from of a specified population from it.
    Returns a pruned tree sequence object.
    """
    tables = ts.dump_tables()
    tables.sites.clear()
    tables.mutations.clear()
    outgroup = ts.get_samples(outgroup_ID)
    """
    Looping through the tree sequence tables and only retaining sites not found in the outgroup.
    """
    for tree in ts.trees():
        for site in tree.sites():
            mut = site.mutations[0]
            pos = site.position
            for leaf in tree.leaves(mut.node):
                if leaf in outgroup:
                    break
            else:
                site_id = tables.sites.add_row(
                position=site.position, ancestral_state=site.ancestral_state)
                tables.mutations.add_row(
                site=site_id, node=mut.node, derived_state=mut.derived_state,parent=mut.parent,time=mut.time)
    return tables.tree_sequence()
