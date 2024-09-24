library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
options(dplyr.print_max = 1e9)


regression=function(df){
	results = data.frame()

    #setting the regression function.

	reg_fun = lm(proportion~date, data = df)

    intercept<-round(coef(reg_fun)[1],3) 
    slope<-round(coef(reg_fun)[2],5)  
    adjust_r2 <-round(as.numeric(summary(reg_fun)[9]),3)

    temp_data_frame = data.frame(SLOPE = slope, INTERCEPT = intercept, R2 = adjust_r2) %>%
    	mutate(label = paste0('intercept = ', INTERCEPT, '  slope = ', SLOPE, '  Adjusted R^2 = ', R2))
    results = rbind(results, temp_data_frame)


    
    results
    }



# -----------------------------------------------------------------------------------------------------------------
# Plot ancient individuals
# -----------------------------------------------------------------------------------------------------------------

exclude = c('Villabruna', 'AfontovaGora3', 'BK1653', 'BKBB7240', 'BKCC7335', 'BKF6620', 'ElMiron', 'OaseNew', 'Ostuni1', 'SalkhitArchAdm', 'Tianyuan', 'Vestonice13', 'Vestonice16', 'Vestonice43', 'LeangPanninge')


data = read.table('Archaic_ancestry.txt',header = T) %>%
	filter(!name %in% exclude, ancient_modern == 'ancient')


model1 = lm(proportion~date, data = data)
summary(model1)



outlier = data %>% filter(date > 20000)


regs = regression(data)

data %>%
	ggplot(aes(x = date, y = proportion)) +
		geom_point() +
		geom_text_repel(data = outlier, aes(x = date , y = proportion, label = name), hjust = 0) +
		geom_text(data = regs, aes(x = 40000, y = 2.5, label = label), hjust = 0, size = 5) +
		geom_smooth(method = 'lm') +
		theme_bw() +
		
		scale_y_continuous('Proportion (%)') +
		scale_x_reverse('Radio Carbon callibrated date')
		



ggsave('3_Archaic_ancestry_overtime_X.pdf', width = 7, height = 5)


# -----------------------------------------------------------------------------------------------------------------
# Plot all
# -----------------------------------------------------------------------------------------------------------------

data = read.table('Archaic_ancestry.txt',header = T) %>%
	filter(!name %in% exclude) %>%
	arrange(region, ancient_modern) 



order = c('EarlyOoA_ancient', 'preLGM-WEurHG_ancient', 'ANS_ancient', 'ANE_ancient', 'Satsurblia_ancient', 'postLGM-WEurHG_ancient', 'SWA_ancient', 'WEur_ancient', 'Yamnaya_ancient', 'postLGM-WEurHG_modern', 'WEur_modern', 'SWA_modern', 'Siberia&Americas_ancient', 'EAS_ancient', 'EAS_modern', 'Siberia&Americas_modern', 'SA_modern', 'Oceania_modern')




skov2023 = data.frame(ID = c('WEur_modern','EAS_modern','Siberia&Americas_modern' ,'SA_modern','Oceania_modern'),
					  proportion = c(0.202 - 0.009, 0.36 - 0.009, 0.237 - 0.027, 0.405 - 0.026 , 0.745 - 0.229))

data %>%
	filter(region != 'Africa') %>%
	mutate(ID = paste(pop, ancient_modern, sep = '_')) %>%
	mutate(ID = factor(ID, order)) %>%
	ggplot() +
		geom_boxplot(aes(x = ID, y = proportion, color = region), alpha = 0.5) +
		geom_jitter(aes(x = ID, y = proportion, color = region), width = 0.10) +
		geom_point(data = skov2023, aes(x = ID, y = proportion), shape = 2, color = 'black') +
		theme_bw() +
		xlab('') + 
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
		



ggsave('4_Archaic_ancestry_X.pdf', width = 7, height = 5)





