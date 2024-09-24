import numpy as np 
from collections import defaultdict
import os
from itertools import combinations


def Get_admixfrog_segments_1000g(name, length_cutoff):

	dict_representation = defaultdict(lambda: np.zeros(300000))
	total_segments, mean_seq_len = 0, []
	
	with open('1000g_hg19_Admixfrogcalls.csv') as data:
		for line in data:
			chrom,start,end,score,target,type,map_len,pos,ID,map_end,pos_end,id_end,seg_len,map_len,pos_len,nscore,n_all_snps,all_n_AFR,all_n_NEA,all_n_DEN,frag_ID,sample = line.strip().split(',')

			if sample == name and int(pos_end) - int(pos) >= length_cutoff:
				pos_start, pos_end = int(pos), int(pos_end)
				start, end = int((pos_start - pos_start % 1000)/1000), int((pos_end - pos_end % 1000)/1000) + 1
				dict_representation[chrom][start:end] = 1
				total_segments += 1
				mean_seq_len.append(pos_end - pos_start)


	return total_segments, np.mean(mean_seq_len), dict_representation



def Get_sstar_segments(name):

	dict_representation = defaultdict(lambda: np.zeros(300000))
	total_segments, mean_seq_len = 0, []

	infile = '1000g_hg19_Sstar.txt'
	with open(infile) as data:
		for line in data:
			ind_name, start, end, length, chrom, region, state, _ = line.strip().split()

			if state == 'neandertal' and ind_name == name:
				pos_start, pos_end = int(start), int(end)
				start, end = int((pos_start - pos_start % 1000)/1000), int((pos_end - pos_end % 1000)/1000) + 1
				dict_representation[chrom][start:end] = 1
				total_segments += 1
				mean_seq_len.append(pos_end - pos_start)


	return total_segments, np.mean(mean_seq_len), dict_representation


def Get_ibdmix_segments(name):

	dict_representation = defaultdict(lambda: np.zeros(300000))
	total_segments, mean_seq_len = 0, []

	infile = '1000g_hg19_IDBMIX.txt'
	with open(infile) as data:
		for line in data:
			if line.startswith('ID'):
				continue

			ind_name, chrom, start, end, LOD, cM = line.strip().split()

			if ind_name == name:
				pos_start, pos_end = int(start), int(end)
				start, end = int((pos_start - pos_start % 1000)/1000), int((pos_end - pos_end % 1000)/1000) + 1
				dict_representation[chrom][start:end] = 1
				total_segments += 1
				mean_seq_len.append(pos_end - pos_start)


	return total_segments, np.mean(mean_seq_len), dict_representation



def Get_CRF_segments_2014(name):

	dict_representation = defaultdict(lambda: np.zeros(300000))
	total_segments, mean_seq_len = 0, []

	infile = '1000g_hg19_CRF_2014.txt'
	with open(infile) as data:
		for line in data:
			ind_name, chrom, start, end = line.strip().split()

			if ind_name == name:
				pos_start, pos_end = int(start), int(end)
				start, end = int((pos_start - pos_start % 1000)/1000), int((pos_end - pos_end % 1000)/1000) + 1
				dict_representation[chrom][start:end] = 1
				total_segments += 1
				mean_seq_len.append(pos_end - pos_start)


	return total_segments, np.mean(mean_seq_len), dict_representation


def Get_hmmix(name):

	dict_representation = defaultdict(lambda: np.zeros(300000))
	total_segments, mean_seq_len = 0, []

	infile = '1000g_hg19_hmmix.txt'
	with open(infile) as data:
		for line in data:
			ind_name, chrom, start, end = line.strip().split()

			if ind_name == name:
				pos_start, pos_end = int(start), int(end)
				start, end = int((pos_start - pos_start % 1000)/1000), int((pos_end - pos_end % 1000)/1000) + 1
				dict_representation[chrom][start:end] = 1
				total_segments += 1
				mean_seq_len.append(pos_end - pos_start)


	return total_segments, np.mean(mean_seq_len), dict_representation



def shared(dict1, dict2, CHROMOSOMES):

	total_dict1 = 0 
	total_dict2 = 0
	total_shared = 0



	for chrom in CHROMOSOMES:
		total_dict1 += np.sum(dict1[chrom])
		total_dict2 += np.sum(dict2[chrom])
		total_shared += np.sum(dict1[chrom] * dict2[chrom])


	return total_shared, total_dict1, total_dict2



CHROMOSOMES = [str(x) for x in range(1,23)]


# 84 CEU individuals
names = ["NA06984", "NA06986", "NA06989", "NA06994", "NA07000", "NA07037", "NA07048", "NA07051", "NA07056", "NA07347", "NA07357", "NA10847", "NA10851", "NA11829", "NA11830", "NA11831", "NA11843", "NA11892", "NA11893", "NA11894", "NA11919", "NA11920", "NA11930", "NA11931", "NA11932", "NA11933", "NA11992", "NA11994", "NA11995", "NA12003", "NA12004", "NA12006", "NA12043", "NA12044", "NA12045", "NA12046", "NA12058", "NA12144", "NA12154", "NA12155", "NA12249", "NA12272", "NA12273", "NA12275", "NA12282", "NA12283", "NA12286", "NA12287", "NA12340", "NA12341", "NA12342", "NA12347", "NA12348", "NA12383", "NA12399", "NA12400", "NA12413", "NA12489", "NA12546", "NA12716", "NA12717", "NA12718", "NA12748", "NA12749", "NA12750", "NA12751", "NA12761", "NA12763", "NA12775", "NA12777", "NA12778", "NA12812", "NA12814", "NA12815", "NA12827", "NA12829", "NA12830", "NA12842", "NA12843", "NA12872", "NA12873", "NA12874", "NA12889", "NA12890"]

summaries = defaultdict(lambda: defaultdict(list))
overlaps = defaultdict(list)

print('name', 'total_shared', 'total_admix', 'admixfrog_segments', 'shared/admixfrog', 'total_other', 'other_segments', 'shared/other', 'other method name', sep = '\t')
for name in names:

	admixfrog_segments, admixfrog_meanlen, admixfrog_dict = Get_admixfrog_segments_1000g(name)


	# overlap with Sstar
	Sstar_segments, Sstar_meanlen, Sstar_dict = Get_sstar_segments(name)
	total_shared, total_admix, total_Sstar = shared(admixfrog_dict, Sstar_dict, CHROMOSOMES)
	print(name, total_shared, total_admix, admixfrog_segments, round(total_shared/total_admix * 100,2), total_Sstar, Sstar_segments, round(total_shared/total_Sstar * 100,2), 'Sstar_2016', sep = '\t')
	summaries['Sstar_2016']['shared'].append(round(total_shared/total_admix * 100,2))

	# overlap with IBDmix
	IBDMIX_segments, IBDMIX_meanlen, IBDMIX_dict = Get_ibdmix_segments(name)
	total_shared, total_admix, total_IBDMIX = shared(admixfrog_dict, IBDMIX_dict, CHROMOSOMES)
	print(name, total_shared, total_admix, admixfrog_segments, round(total_shared/total_admix * 100,2), total_IBDMIX, IBDMIX_segments, round(total_shared/total_IBDMIX * 100,2), 'IBDMIX_2020', sep = '\t')
	summaries['IBDMIX_2020']['shared'].append(round(total_shared/total_admix * 100,2))

	# overlap with CRF 2014
	CRF_segments, CRF_meanlen, CRF_dict = Get_CRF_segments_2014(name)
	total_shared, total_admix, total_CRF = shared(admixfrog_dict, CRF_dict, CHROMOSOMES)
	print(name, total_shared, total_admix, admixfrog_segments, round(total_shared/total_admix * 100,2), total_CRF, CRF_segments, round(total_shared/total_CRF * 100,2), 'CRF_2014', sep = '\t')
	summaries['CRF_2014']['shared'].append(round(total_shared/total_admix * 100,2))

	# overlap with hmmix
	hmmix_segments, hmmix_meanlen, hmmix_dict = Get_hmmix(name)
	total_shared, total_admix, total_hmmix = shared(admixfrog_dict, hmmix_dict, CHROMOSOMES)
	print(name, total_shared, total_admix, admixfrog_segments, round(total_shared/total_admix * 100,2), total_hmmix, hmmix_segments, round(total_shared/total_hmmix * 100,2), 'hmmix_2018', sep = '\t')
	summaries['hmmix_2018']['shared'].append(round(total_shared/total_admix * 100,2))

	print()

	summaries['admixfrog']['meanlen'].append(admixfrog_meanlen)
	summaries['admixfrog']['segments'].append(admixfrog_segments)
	summaries['admixfrog']['total'].append(total_admix)
	summaries['admixfrog']['shared'].append(1)

	summaries['CRF_2014']['meanlen'].append(CRF_meanlen)
	summaries['CRF_2014']['segments'].append(CRF_segments)
	summaries['CRF_2014']['total'].append(total_CRF)
	
	summaries['hmmix_2018']['meanlen'].append(hmmix_meanlen)
	summaries['hmmix_2018']['segments'].append(hmmix_segments)
	summaries['hmmix_2018']['total'].append(total_hmmix)

	summaries['IBDMIX_2020']['meanlen'].append(IBDMIX_meanlen)
	summaries['IBDMIX_2020']['segments'].append(IBDMIX_segments)
	summaries['IBDMIX_2020']['total'].append(total_IBDMIX)
	
	summaries['Sstar_2016']['meanlen'].append(Sstar_meanlen)
	summaries['Sstar_2016']['segments'].append(Sstar_segments)
	summaries['Sstar_2016']['total'].append(total_Sstar)


	all_dicts = {'admixfrog': admixfrog_dict,
				'hmmix': hmmix_dict,
				'CRF': CRF_dict,
				'Sstar': Sstar_dict,
				'IBDMIX': IBDMIX_dict}

	for (a, b) in combinations(all_dicts.keys(), 2):
		total_shared, total_a, total_b = shared(all_dicts[a], all_dicts[b], CHROMOSOMES)
		#print(a, b, total_shared, total_a, total_b)

		overlaps[f'{a}|{b}'].append(total_shared/total_a)
		overlaps[f'{b}|{a}'].append(total_shared/total_b)




for key, values in overlaps.items():
	a, b = key.split('|')
	print(a, b, np.mean(values), sep = '\t')
print()



for method in summaries:
	for stat in ['meanlen', 'segments', 'total', 'shared']:
		print(method, stat, np.mean(summaries[method][stat]), sep = '\t')
	print()










