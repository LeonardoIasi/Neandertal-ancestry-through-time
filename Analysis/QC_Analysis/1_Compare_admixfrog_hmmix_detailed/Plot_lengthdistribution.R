library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)



# --------------------------------------------------------------------------------------------------------------------
# Length distribution and example plot


# load all simulated files
data_lengthdist = read.table('hmmix_admixfrog_detailed.txt', header = T) %>%
	mutate(ycoord = ifelse(method == 'hmmix', 0, 1))


chromosomes = data.frame(chrom = c(1:22),
						starts = rep(0,22),
						ends = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566)) %>%
	mutate(chrom = factor(chrom, 1:22) ) %>%
	mutate(labels = paste0('chr', chrom, ' '))


gaps = read.table('gaps.txt', header = T) %>%
	filter(!chrom %in% c('X','Y'), size > 1e6) %>%
	mutate(chrom = factor(chrom, 1:22) )

data_lengthdist %>%
	#filter(name == 'B_French-3') %>%
	filter(name == 'S_Aleut-2') %>%
	mutate(chrom = factor(chrom, 1:22)) %>%

	ggplot() +


	# chromosomes positions
	geom_rect(data = chromosomes, aes(xmin = starts/1e6, xmax = ends/1e6, ymin = -0.02, ymax = 2.02), fill = NA, color = 'black') +


	geom_text(data = chromosomes, aes(x = 0, y = 1, label = labels), hjust = 1) +
	geom_rect(data = gaps, aes(xmin = start/1e6, xmax = end/1e6, ymin = 0.02, ymax = 1.99), fill = 'darkgrey') +


	geom_rect(aes(xmin = start/1e6, xmax = end/1e6, ymin = ycoord, ymax = ycoord + 1, fill = method)) +

	geom_rect(data = chromosomes, aes(xmin = starts/1e6, xmax = ends/1e6, ymin = -0.02, ymax = 2.02), fill = NA, color = 'black') +


	facet_grid(chrom~., switch = 'y') +
	theme_bw() +
	xlab('Genomic position (Mb)') + 
	scale_y_continuous('', expand = c(0,0)) +
	labs(x = 'genomic position (Mb)') +
		  theme(strip.background = element_blank(),
	        axis.title.y=element_blank(),
	        axis.text.y=element_blank(),
	        axis.ticks.y=element_blank(),
	        panel.grid.major = element_blank(),
	      panel.grid.minor = element_blank(),
	      panel.border = element_blank(),
	      panel.background = element_blank(),
	      strip.text.y = element_blank(), #element_text(angle = 180)
	      legend.position = 'bottom',
	      panel.spacing.y=unit(0, "lines"))


ggsave('0_Chromosome_view.pdf', width = 15, height = 7)

name_order = c('B_Karitiana-3', 'B_Mixe-1', 'S_Chane-1', 'S_Pima-1', 'S_Pima-2', 'S_Karitiana-1', 'S_Karitiana-2', 'S_Mayan-1', 'S_Mayan-2', 'S_Mixe-2', 'S_Altaian-1', 'S_Chukchi-1', 'S_Eskimo_Chaplin-1', 'S_Itelman-1', 'S_Aleut-1', 'S_Aleut-2', 'S_Eskimo_Naukan-1', 'S_Eskimo_Naukan-2', 'S_Eskimo_Sireniki-1', 'S_Kyrgyz-1', 'B_Dai-4', 'B_Han-3', 'S_Atayal-1', 'S_Daur-2', 'S_Japanese-1', 'S_Ami-1', 'S_Ami-2', 'S_Burmese-1', 'S_Burmese-2', 'S_Cambodian-1', 'B_Australian-3', 'B_Australian-4', 'B_Papuan-15', 'S_Hawaiian-1', 'S_Maori-1', 'S_Papuan-3', 'S_Bougainville-1', 'S_Bougainville-2', 'S_Dusun-1', 'S_Dusun-2', 'S_Igorot-1', 'S_Khonda_Dora-1', 'S_Balochi-1', 'S_Balochi-2', 'S_Punjabi-4', 'S_Bengali-1', 'S_Bengali-2', 'S_Brahmin-1', 'S_Brahmin-2', 'S_Brahui-1', 'S_Brahui-2', 'S_Burusho-1', 'B_Crete-1', 'B_Crete-2', 'B_French-3', 'B_Sardinian-3', 'S_Albanian-1', 'S_Chechen-1', 'S_Czech-2', 'S_Norwegian-1', 'S_Polish-1', 'S_Samaritan-1', 'S_Armenian-1')

data_lengthdist %>%
	#filter(segmentlength < 2e6) %>%
	mutate(name = factor(name, name_order)) %>%
	ggplot() +
	geom_histogram(aes(x = segmentlength), bins = 50) +
	facet_grid(method~.) +
	scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
	scale_y_log10() + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	annotation_logticks() 


ggsave('1_Lengthdistribution.pdf', width = 15, height = 7)


data_lengthdist %>%
	group_by(name, region, method) %>%
	summarize(mean_length = mean(segmentlength), total_segments = n()) %>%
	ungroup() %>%
	gather(key, value, mean_length:total_segments) %>%
	mutate(name = factor(name, name_order)) %>%
	ggplot() +
	geom_point(aes(x = name, y = value, color = method)) +
	facet_grid(key~., scale = 'free_y') +
	theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


ggsave('2_Lengthdistribution_summary.pdf', width = 15, height = 7)


data_lengthdist %>%
	group_by(method, name) %>%
	summarize(mean_length = mean(segmentlength), total_segments = n(), total_sequence = sum(segmentlength)) %>%
	group_by(method) %>%
	summarize(av_mean_length = mean(mean_length), av_total_seq = mean(total_sequence), av_segments = mean(total_segments))

	


# --------------------------------------------------------------------------------------------------------------------
# Detailed overlap analysis - remove long segments so it is not affected by admixfrog spanning centromeres


# Plot detailed
data = data_lengthdist %>%
	filter(segmentlength < 1e6)




data %>%
	mutate(rounded_overlap = overlap_ratio - overlap_ratio%%0.05) %>%
	group_by(method, name, rounded_overlap) %>%
	summarize(total_overlap = sum(segmentlength)) %>%
	ungroup() %>%
	ggplot() +
	geom_jitter(aes(x = rounded_overlap, y = total_overlap/1e6 ), width = 0.01) +
	geom_boxplot(aes(x = rounded_overlap, y = total_overlap/1e6, group = rounded_overlap )) +
	facet_grid(method~.) +
	theme_bw() +
	ylab('Total_sequence (Mb)') +
	xlab('Overlap')


ggsave('3_Overlapratio_sequence.pdf')


data %>%
	mutate(rounded_overlap = overlap_ratio - overlap_ratio%%0.05) %>%
	group_by(method, rounded_overlap) %>%
	summarize(mean_length = mean(segmentlength)) %>%
	ungroup() %>%
	ggplot() +
	geom_point(aes(x = rounded_overlap, y = mean_length)) +
	facet_grid(method~.) +
	theme_bw()


ggsave('4_Overlapratio_length.pdf')




print('hmmix unique')

data %>%
	filter(method == 'hmmix') %>%
	rowwise() %>%
	mutate(arch_snps = max(altai, vindija, chagyrskaya)) %>%
	mutate(rounded_overlap = ifelse(overlap_ratio > 0, 'shared', 'not_shared') ) %>%
	group_by(rounded_overlap, name) %>%
	summarize(total_length = sum(segmentlength), n_segments = n(), shared = mean(arch_snps/allsnps)) %>%
	group_by(rounded_overlap) %>%
	summarize(av_seg_len = mean(total_length), av_match_rate = mean(shared), av_segments = mean(n_segments)) 




data %>%
	filter(method == 'hmmix') %>%
	rowwise() %>%
	mutate(arch_snps = max(altai, vindija, chagyrskaya)) %>%
	mutate(snp_sum = arch_snps/allsnps) %>%
	mutate(overlaps = ifelse(overlap_ratio > 0, 'shared', 'not_shared') ) %>%
	

	mutate(rounded_snps_sim = snp_sum - snp_sum%%0.05) %>%
	group_by(rounded_snps_sim, overlaps, name) %>%
	summarize(total_length = sum(segmentlength), n_segments = n()) %>%
	group_by(rounded_snps_sim, overlaps) %>%
	summarize(av_total_length = mean(total_length)/1e6, av_n_segments = mean(n_segments)) %>%
	ungroup() %>%
	gather(key, value, av_total_length:av_n_segments) %>%

	mutate(overlaps = factor(overlaps, c('shared', 'not_shared'))) %>%
	ggplot() +
	geom_bar(aes(x = rounded_snps_sim, y = value, fill = overlaps), stat = 'identity', position = 'stack') +
	facet_wrap(~key, scale = 'free_y', ncol = 2) +
	xlab('Match rate to any Neanderthal') +
	ylab('count/sequence') +
	theme_bw()


ggsave('5.1_Overlap_archaic.pdf', width = 15, height = 7)


print('segments length')

data %>%
	filter(segmentlength < 1e6) %>%
	mutate(overlaps = ifelse(overlap_ratio > 0, 'shared', 'not_shared') ) %>%
	group_by(overlaps, name) %>%
	summarize(mean_length = mean(segmentlength)) %>%
	group_by(overlaps) %>%
	summarize(av_seg_len = mean(mean_length)) 



data %>%
	filter(segmentlength < 500000) %>%
	mutate(overlaps = ifelse(overlap_ratio > 0, 'shared', 'not_shared') ) %>%
	
	mutate(rounded_snps_sim = segmentlength - segmentlength%%10000) %>%
	group_by(method, rounded_snps_sim, overlaps, name) %>%
	summarize(total_length = sum(segmentlength), n_segments = n()) %>%
	group_by(method, rounded_snps_sim, overlaps) %>%
	summarize(av_total_length = mean(total_length)/1e6, av_n_segments = mean(n_segments)) %>%
	ungroup() %>%
	gather(key, value, av_total_length:av_n_segments) %>%

	mutate(overlaps = factor(overlaps, c('shared', 'not_shared'))) %>%
	ggplot() +
	geom_bar(aes(x = rounded_snps_sim, y = value, fill = overlaps), stat = 'identity', position = 'stack') +
	facet_wrap(method~key, scale = 'free_y', ncol = 2) +
	xlab('Segment length') +
	ylab('count/sequence') +
	theme_bw()


ggsave('5.2_Overlap_length.pdf', width = 15, height = 7)

old_way = data %>%
	filter(segmentlength < 250000, segmentlength > 10000) %>%
	mutate(overlaps = ifelse(overlap_ratio > 0, 'shared', 'not_shared') ) %>%
	

	mutate(rounded_snps_sim = segmentlength - segmentlength%%5000) %>%
	group_by(method, rounded_snps_sim, name) %>%
	mutate(total_sequence = sum(segmentlength), total_segments = n()) %>%
	group_by(method, rounded_snps_sim, overlaps, name) %>%
	summarize(shared_length = sum(segmentlength)/mean(total_sequence), shared_segments = n()/mean(total_segments)) %>%
	group_by(method, rounded_snps_sim, overlaps) %>%
	filter(overlaps == 'shared') %>%
	summarize(av_total_length = mean(shared_segments), conf = 1.96 * sd(shared_segments)/sqrt(n()) ) %>%
	mutate(MIN = av_total_length - conf, MAX = av_total_length + conf) %>%
	ungroup() %>%
	mutate(facet = 'catagorical shared/not shared') %>%
	select(-overlaps)




new_data = data %>%
	filter(segmentlength < 250000, segmentlength > 10000) %>%
	mutate(overlaps = overlap_ratio) %>%
	

	mutate(rounded_snps_sim = segmentlength - segmentlength%%5000) %>%
	group_by(method, rounded_snps_sim, name) %>%
	mutate(shared_segments = sum(segmentlength * overlap_ratio) / sum(segmentlength) ) %>%
	group_by(method, rounded_snps_sim) %>%
	summarize(av_total_length = mean(shared_segments), conf = 1.96 * sd(shared_segments)/sqrt(n()) ) %>%
	mutate(MIN = av_total_length - conf, MAX = av_total_length + conf) %>%
	ungroup() %>%
	mutate(facet = 'shared seg/total sequence')

rbind(new_data, old_way) %>%
	ggplot(aes(x = rounded_snps_sim, y = av_total_length * 100, color = facet)) +
	geom_ribbon(aes(ymin = MIN * 100, ymax = MAX * 100), fill = "grey70") +
	geom_line() +
	facet_grid(method~.) +
	xlab('Segment length') +
	scale_y_continuous('% Shared') +
	coord_cartesian(ylim = c(0,100)) +
	theme_bw()


ggsave('5.3_Overlap_length_percent.pdf', width = 12, height = 4)




# Why are hmmix not finding those in admixfrog? They are shorter and have fewer archaic snps
print('Admixfrog unique')
admix_frog_data = read.table('hmmix_admixfrog_detailed_SNPs.txt', header = T)

admix_frog_data %>%
	mutate(snp_sum = derived_snps/segmentlength * 1000) %>%
	mutate(rounded_overlap = ifelse(overlap_ratio > 0, 'shared', 'not_shared') ) %>%
	group_by(rounded_overlap, name) %>%
	summarize(average_derived_snps_name = mean(snp_sum), average_length_name = mean(segmentlength), n_segments_name = n(), total_length = sum(segmentlength)) %>%
	group_by(rounded_overlap) %>%
	summarize(average_derived_snps = mean(average_derived_snps_name), average_length = mean(average_length_name), n_segments = mean(n_segments_name), av_seq = mean(total_length)) %>%
	ungroup() 


admix_frog_data %>%
	mutate(snp_sum = derived_snps/segmentlength * 1000) %>%
	filter(snp_sum < 1) %>%
	mutate(overlaps = ifelse(overlap_ratio > 0, 'shared', 'not_shared') ) %>%
	mutate(overlaps = factor(overlaps, c('shared', 'not_shared'))) %>%
	ggplot() +
	geom_histogram(aes(x = snp_sum, fill = overlaps)) +
	geom_vline(aes(xintercept = 0.35 * 0.75), linetype = 'dashed') +
	#facet_grid(overlaps~.) +
	theme_bw() +
	xlab('Derived snps per 1kb')

ggsave('6_derived_density.pdf')



admix_frog_data %>%
	mutate(snp_sum = derived_snps/segmentlength * 1000) %>%
	filter(snp_sum < 1) %>%
	mutate(overlaps = ifelse(overlap_ratio > 0, 'shared', 'not_shared') ) %>%
	mutate(overlaps = factor(overlaps, c('shared', 'not_shared'))) %>%
	ggplot() +
	geom_histogram(aes(x = snp_sum, fill = overlaps)) +
	geom_vline(aes(xintercept = 0.35 * 0.75), linetype = 'dashed') +
	#facet_grid(overlaps~.) +
	theme_bw() +
	xlab('Derived snps per 1kb')

ggsave('7_derived_density.pdf')




admix_frog_data %>%
	rowwise() %>%
	mutate(snp_sum = derived_snps/segmentlength * 1000) %>%
	filter(snp_sum < 1) %>%
	mutate(overlaps = ifelse(overlap_ratio > 0, 'shared', 'not_shared') ) %>%
	mutate(rounded_snps_sim = snp_sum - snp_sum%%0.05) %>%
	group_by(rounded_snps_sim, overlaps, name) %>%
	summarize(total_length = sum(segmentlength), n_segments = n()) %>%
	group_by(rounded_snps_sim, overlaps) %>%
	summarize(av_total_length = mean(total_length)/1e6, av_n_segments = mean(n_segments)) %>%
	ungroup() %>%
	gather(key, value, av_total_length:av_n_segments) %>%
	mutate(overlaps = factor(overlaps, c('shared', 'not_shared'))) %>%
	ggplot() +
	geom_bar(aes(x = rounded_snps_sim, y = value, fill = overlaps), stat = 'identity', position = 'stack') +
	geom_vline(aes(xintercept = 0.35 * 0.75), linetype = 'dashed') +
	facet_wrap(~key, scale = 'free_y', ncol = 2) +
	xlab('SNP density') +
	ylab('count/sequence') +
	theme_bw()


ggsave('7.1_Overlap_archaic.pdf', width = 15, height = 7)