library(ggplot2)
library(dplyr)
options(dplyr.print_max = 1e9)

# Plot overlap
data = read.table('hmmix_admixfrog.txt', header = T, sep = '\t') %>%
	filter(!name %in% c('S_Eskimo_Chaplin-1', 'S_Eskimo_Naukan-1', 'S_Eskimo_Naukan-2', 'S_Eskimo_Sireniki-1'))



data %>%
	mutate(percent_found = shared/total * 100) %>%
	group_by(archaic, method) %>%
	#filter(method == 'hmmix') %>%
	summarize(average = mean(percent_found, na.rm = TRUE))



data %>%
	mutate(percent_found = shared/total * 100) %>%
	filter(archaic != 'BOTH') %>%
	group_by(archaic, region, method) %>%
	summarize(method_found = mean(total), average = mean(percent_found, na.rm = TRUE)) %>%
	arrange(region, archaic, method)




to_keep = data %>%
	filter(archaic == 'BOTH', method == 'hmmix') %>%
	group_by(region) %>%
	sample_n(10) %>%
	select(name)

data %>%
	filter(archaic != 'BOTH', name %in% to_keep$name) %>%
	ggplot() +
	geom_bar(aes(x = name, y = total/1000, fill = method), position = 'dodge', stat = 'identity', alpha = 0.8) +
	geom_bar(aes(x = name, y = shared/1000, fill = method),  color = 'black', position = 'dodge', stat = 'identity') +
	facet_wrap(archaic~region, scale = 'free', nrow = 2) +
	scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
	theme_bw() +
	theme(panel.spacing.x=unit(0.0, "lines"), 
        panel.spacing=unit(0,"lines"),
        
        axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	ylab('Called sequence (Mb)') +
	ggtitle('Total amount of called sequence (Black bars indicates sharing)')

ggsave('8_Admixfrog_compare_to_hmmix.pdf', width = 12)


