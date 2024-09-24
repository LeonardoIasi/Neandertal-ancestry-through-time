library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
options(dplyr.print_max = 1e9)

exclude = c('Villabruna', 'AfontovaGora3', 'BK1653', 'BKBB7240', 'BKCC7335', 'BKF6620', 'ElMiron', 'OaseNew', 'Ostuni1', 'SalkhitArchAdm', 'Tianyuan', 'Vestonice13', 'Vestonice16', 'Vestonice43', 'LeangPanninge')


BORDER = 250000


regions = read.table('X_sweeps.txt', col.names = c('START', 'END'))


not_callable = read.table('X_ArchaicAdmixtureX_maske_file.bed', col.names = c('chrom','start', 'end')) %>%
	mutate(len = end - start) %>%
	filter(len > 20000)




Segments = read.table('segments_X.txt', header = T) %>%
	filter(!name %in% exclude) %>%
	mutate(facet_plot = ifelse(ancient_modern == 'modern', as.character(region), 'ancient') ) %>%
	mutate(facet_plot = ifelse(date > 44000, 'EarlyOoA', facet_plot) ) 


row_numbers = Segments %>%
	distinct(name, facet_plot, region, ancient_modern, date) %>%
	arrange(date, facet_plot) %>%
	mutate(rownumber = 1:n()) %>%
	ungroup()




uncalled = data.frame(facet_plot = 'uncalled', yvalue = -10, label_value = -5)


break_points = row_numbers %>%
	group_by(facet_plot) %>%
	summarize(yvalue = min(rownumber) - 0.3, label_value = mean(rownumber) - 0.3) %>%
	ungroup() %>%
	rbind(uncalled)


# Entire genome
recomb_map = read.table('maps_chr.X', header = T)


rate_diffs = recomb_map %>%
	select(Physical_Pos, AA_Map) %>%
	mutate(previous_pos = lag(Physical_Pos), previous_pos_rec = lag(AA_Map)) %>%
	mutate(pos_change = (Physical_Pos - previous_pos)/1e6, rec_change = AA_Map - previous_pos_rec ) %>%
	mutate(rate = rec_change/pos_change) %>%
	mutate(rounded_pos = Physical_Pos - Physical_Pos %% 20000) %>%
	group_by(rounded_pos) %>%
	summarize(mean_rate = mean(rate, na.rm = TRUE))



rec_vs_freq = read.table('Rec_rate_vs_freq_ancient_modern.txt', header = T, sep = '\t')


p1 = rec_vs_freq %>%
	ggplot() +
	geom_line(aes(x = start + 10000, y = freq_modern), color = "#56B4E9") +
	geom_line(aes(x = start + 10000, y = -freq_ancient), color = "#E69F00") +
	theme_bw() +
	xlab('') +
	theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
	#ylab() 
	scale_y_continuous('Frequency (%)', breaks = c(-10, -5, 0, 5, 10), labels =c(10, 5, 0, 5, 10)) #+coord_cartesian(ylim = c(0, 100))


p2 = rate_diffs %>%
	ggplot() +
	geom_line(aes(x = rounded_pos, y = mean_rate)) +
	theme_bw() +
	xlab('') +
	theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
	#ylab() 
	scale_y_continuous('Recombination rate\n(cM/Mb)') +
	coord_cartesian(ylim = c(0, 100))



p3 = Segments %>%
	inner_join(row_numbers, by = c('name', 'region', 'ancient_modern')) %>%
	ggplot() +

		# As positions
		geom_rect(data = regions, aes(xmin = START/1e6, xmax = END/1e6, ymin = 0.7, ymax = Inf ), fill = 'lightgrey') +
		geom_rect(data = not_callable, aes(xmin = start/1e6, xmax = end/1e6, ymin = -10, ymax = 0.7 ), color = 'darkred', fill = "darkred") +
		geom_rect(aes(xmin = start/1e6, xmax = end/1e6, ymin = rownumber - 0.1, ymax = rownumber + 0.1, fill = ancient_modern, color = ancient_modern), linewidth = 1) +
		

		theme_bw() +
		scale_y_continuous('', breaks = break_points$label_value, labels = break_points$facet_plot, expand = c(0,0)) +
		xlab('Genomic position (Mb)') +

		scale_fill_manual('Region',values = c("#E69F00", "#56B4E9")) +
		scale_color_manual('Region',values = c("#E69F00", "#56B4E9")) +
		
		geom_hline(data = break_points, aes(yintercept = yvalue)) +

		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
			  panel.background = element_blank(), axis.line = element_line(colour = "black"), 
			  panel.spacing = unit(0,'lines'),
			  strip.text.y.left = element_text(angle = 0),
			  legend.position = 'bottom')

plot_grid(p2, p1, p3, align = "hv", ncol = 1, rel_heights = c(0.15, 0.15, 1))
ggsave('1_Entire_X_sweeps.pdf', width = 12, height = 18)








