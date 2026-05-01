
#############################
## LIBRARIES

library(ggplot2)
library(ggVennDiagram)
library(here)
library(xtable)

# also requires:
# cowplot, reshape2, stringr

######################
## READ DATA

# methods and nicer string conversion
methods <- c('cisTWAS', 'cisPWAS', 'cistransTWAS', 'cistransPWAS')
method_titles <- setNames(c('cis-TWAS', 'cis-PWAS', 'cis/trans-TWAS', 'cis/trans-PWAS'), methods)

# column names
dat_cols <- c('Pathway_Source', 'Pathway_Name', 
	'pvalue', 'qvalue_FDR_BH', 'qvalue_Bonferroni', 
	'Type', 'Category', 'method')

# read in data
dat <- data.frame()

for (i in 1:4) {
	# read data
	dat_raw <- read.table(
		here('data', 
			paste0(methods[i], '_pathDIP_Enrichment_analysis_for_Panther_Pathway.txt')), 
		header = TRUE, 
		sep = '\t', 
		skip = 2, # pathDIP output adds comment lines to first 2 rows
		check.names = FALSE)
	dat_raw$method <- methods[i]

	# set column names
	colnames(dat_raw) <- dat_cols

	dat <- rbind(dat, dat_raw)
}

# remove row names
row.names(dat) <- NULL
dat$neglogFDR <- -log10(dat$qvalue_FDR_BH)

# reduce length of strings by removing " pathway" if it occurs at end of pathway name
dat$Pathway_Name_old <- dat$Pathway_Name
dat$Pathway_Name <- gsub(' pathway$', '', dat$Pathway_Name_old)


######################
## PATH PLOTS

## function to make path plot
path_plot <- function(dat, title, twidth = 43, text_size = 16){
	dat <- dat[order(-dat$qvalue_FDR_BH), ]

	dat$path <- factor(dat$Pathway_Name, 
		levels = unique(dat$Pathway_Name), 
		labels = stringr::str_wrap(unique(dat$Pathway_Name), twidth))

	# assign colors
	pal <- c('#3A948E', '#36807A', '#2f615d', '#44B5AD')
	plot_colors <- rep(pal, 
		(nrow(dat) %/% 4) + 
			as.integer( (nrow(dat) %% 4) != 0 ) )[1:nrow(dat)]

	# create ggplot
	p <- ggplot(data = dat, 
			aes(x = path, y = neglogFDR, 
				color = path, fill = path)
		) + 
		geom_col(alpha = 0.95) +
		labs(
			title = title, x = '', 
			y = bquote(bold(-"log"["10"]("FDR")))
		) +
		scale_color_manual(values = plot_colors) +
		scale_fill_manual(values = plot_colors) +
		theme_bw() +
		guides(color = 'none', fill = 'none') +
		theme(
			text = element_text(
				size = text_size+1, 
				face = 'bold'),
			axis.text.y = element_text(
				size = text_size+1, 
				color = 'black', face = 'bold'),
			axis.text.x = element_text(
				size = text_size+1, 
				color = 'black', face = 'bold'),
			axis.title.x = element_text(
				size = text_size+3, 
				color = 'black', face = 'bold'),
			panel.grid.major.y = element_blank(),
			panel.grid.minor.y = element_blank()
		) +
		scale_y_continuous(
			breaks = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5), 
			limits = c(0,2.6)
		) +
		coord_flip()

	return(p)
}

# function to get plot for one method
method_path_plot <- function(method){
	p <- path_plot(dat[dat[['method']] == method, ], method_titles[[method]])
	return(p)
}

paths_plot_list <- lapply(methods, method_path_plot)

# grid of paths plots, one for each method
paths_plot <- cowplot::plot_grid(
	plotlist = paths_plot_list, 
	labels = paste0(LETTERS[1:4],'.'),
	label_size = 28)

ggsave(here('figs', 'paths_plot.pdf'), paths_plot, width = 24, height = 26)


#################
# VENN DIAGRAM OVERLAPPED PATHS

# factor method
dat$method <- factor(dat$method, 
	levels=c('cistransTWAS', 'cisTWAS', 'cisPWAS', 'cistransPWAS'),
	labels=c('cis/trans-TWAS', 'cis-TWAS', 
		'cis-PWAS', 'cis/trans-PWAS'))

# get list for venn diagram
vennlist <- as.list(by(dat$Pathway_Name, dat$method, function(x) x))

venn_plot <- ggVennDiagram(vennlist, label_alpha=0, label='count') + 
	theme(legend.pos='None') + 
	scale_fill_gradient(low='white', high='#36908a') +
	scale_x_continuous(expand = expansion(mult=0.15)) # expand so long labels fit

ggsave(here('figs', 'paths_venn_report.pdf'), venn_plot)


#################
# TABLE OVERLAPPED PATHS

## overlapping paths
path_dat <- reshape2::dcast(dat, Pathway_Name ~ method, length, value.var='Pathway_Name')

# reorder cols
path_dat <- path_dat[, c('Pathway_Name', 'cis-TWAS', 'cis/trans-TWAS',  
		'cis-PWAS', 'cis/trans-PWAS')]

# total number of methods path id'd with
path_dat$total <- rowSums(path_dat[, c('cis-TWAS', 'cis/trans-TWAS',  
		'cis-PWAS', 'cis/trans-PWAS')])

# get total number for each xwas type
path_dat$twas_total <- rowSums(path_dat[, c('cis-TWAS', 'cis/trans-TWAS')])
path_dat$pwas_total <- rowSums(path_dat[, c('cis-PWAS', 'cis/trans-PWAS')])

# sort
path_dat <- path_dat[ with(path_dat, order(-total, -twas_total, -pwas_total, -`cis-TWAS`, -`cis/trans-TWAS`, 
		-`cis-PWAS`, -`cis/trans-PWAS`)), ]

# remove rownames
row.names(path_dat) <- NULL

# table for report
all_overlap_xtable <- xtable(
	path_dat[path_dat$total == 4, 'Pathway_Name', drop = FALSE], 
	display = rep('s', 2)
	)

caption(all_overlap_xtable) <- 'Significantly Enriched Pathways for all xWAS'

label(all_overlap_xtable) <- 'tab:allenriched'

print(all_overlap_xtable, 
	include.colnames = FALSE,
	include.rownames = FALSE, 
	booktabs = TRUE, 
	table.placement = 'H', 
	sanitize.text.function = function(x){x},
	caption.placement = 'top',
	append = FALSE,
	file = here('tables', 'overlap_paths.tex'))




