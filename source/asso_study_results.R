
#############################
## LIBRARIES

# also requires: here

# all plots
library(ggplot2)

# Manhattan plot
library(ggrepel)
library(ggnewscale)

# Venn diagram
library(ggVennDiagram)

# here function
library(here)

# dataframe manipulation
library(reshape2)

# LaTeX tables
library(xtable)


#############################
## SOURCE

# Manhattan plot script
source(here('source', 'manhattan_plot.R'))

#############################
## ANNOTATION

## function to get annotation dataframe
get_annot_data <- function(filename){
	annot <- read.csv(here('data', filename),
		sep='\t', header=TRUE	,
		colClasses=c('character', rep('integer', 2), 'character')
	)

	# trim whitespace frome gene name
	annot$GeneName <- trimws(as.character(annot$GeneName))

	# should only have one row per GeneName
	annot <- annot[!duplicated(annot$GeneName), ]

	return(annot)
}

## get annotation data for PWAS and TWAS
annot_pwas <- get_annot_data('ProteinAbundance_ROSMAP_annot.txt')
annot_twas <- get_annot_data('TranscriptAbundance_ROSMAP_annot.txt')

## function to annotate data
annotate_dat <- function(data, twas){
	if (twas == 'PWAS') {
		annot  <- annot_pwas
	} else 	if (twas == 'TWAS') {
		annot  <- annot_twas
	}

	# trim GeneName whitespace
	data$GeneName <- trimws(as.character(data$GeneName))

	# annotate results
	data <- merge(annot, data, by='GeneName', all.x=FALSE, all.y=TRUE)

	# reorder by chrom, genestart
	data <- data[order(data$CHROM, data$GeneStart),]

	return(data)
}


#############################
## ADJUST PVAL FUNCTION

# adjust for genomic control factor
chisq_0.5 <- qchisq(1/2, df=1)

GC_adjust <- function(pval) {
	pval <- as.numeric(pval)
	chisq_stats <- qchisq(pval, df=1, lower.tail=FALSE)

	# lambda
	GC <- median(chisq_stats, na.rm=TRUE) / chisq_0.5
	print(GC)

	adj_chisq <- chisq_stats / GC
	adj_p <- pchisq(adj_chisq, df=1, lower.tail=FALSE)

	return(cbind(chisq_stats, adj_chisq, adj_p, GC))
}


#############################
## RESULTS FILES

# directory
dat_cols <- c('GeneName', 'Z', 'PVAL', 'N_with_weight', 
	'Ncis_with_weight', 'N_tested', 'Ncis_tested', 
	'sum_cpp_cis', 'sum_cpp_trans')

# cis-PWAS
c_pwas <- read.csv(here('data', 'cis_PWAS_results.txt'),
	sep='\t', header=TRUE)
colnames(c_pwas) <- dat_cols
c_pwas$snp_type <- 'cis'
c_pwas$method_type <- 'PWAS'
c_pwas$method <- 'cis-PWAS'
c_pwas[, c('chisq', 'adj_chisq', 'adj_PVAL', 'lambda')] <- GC_adjust(c_pwas$PVAL)
c_pwas$zero_cis <- FALSE

# cis/trans-PWAS
ct_pwas <- read.csv(here('data', 'cistrans_PWAS_results.txt'),
	sep='\t', header=TRUE)
colnames(ct_pwas) <- c(dat_cols, 'zero_cis')
ct_pwas$snp_type <- 'cis/trans'
ct_pwas$method_type <- 'PWAS'
ct_pwas$method <- 'cis/trans-PWAS'
ct_pwas[, c('chisq', 'adj_chisq', 'adj_PVAL', 'lambda')] <- GC_adjust(ct_pwas$PVAL)

# cis-TWAS
c_twas <- read.csv(here('data', 'cis_TWAS_results.txt'),
	sep='\t', header=TRUE)
colnames(c_twas)[1:3] <- dat_cols
c_twas$snp_type <- 'cis'
c_twas$method_type <- 'TWAS'
c_twas$method <- 'cis-TWAS'
c_twas[, c('chisq', 'adj_chisq', 'adj_PVAL', 'lambda')] <- GC_adjust(c_twas$PVAL)
c_twas$zero_cis <- FALSE

# cis/trans-TWAS
ct_twas <- read.csv(here('data', 'cistrans_TWAS_results.txt'),
	sep='\t', header=TRUE)
colnames(ct_twas)[1:3] <- dat_cols
ct_twas$snp_type <- 'cis/trans'
ct_twas$method_type <- 'TWAS'
ct_twas$method <- 'cis/trans-TWAS'
ct_twas[, c('chisq', 'adj_chisq', 'adj_PVAL', 'lambda')] <- GC_adjust(ct_twas$PVAL)
ct_twas$zero_cis <- ct_twas$Ncis_tested == 0


#############################
## COMBINE RESULTS: combine dat for all methods

# combine cis-PWAS, cis/trans-PWAS
dat_pwas <- rbind(c_pwas, ct_pwas)
dat_pwas <- annotate_dat(dat_pwas, 'PWAS')

# combine dat with cis-TWAS, cis/trans-TWAS
dat_twas <- rbind(c_twas, ct_twas)
dat_twas <- annotate_dat(dat_twas, 'TWAS')

cols_keep <- c('CHROM', 'GeneStart', 'GeneEnd', 'GeneName',
	'Z', 'chisq', 'PVAL', 'adj_chisq', 'adj_PVAL', 'lambda',
	'N_with_weight', 'Ncis_with_weight', 'N_tested', 'Ncis_tested', 
	'sum_cpp_cis', 'sum_cpp_trans', 'zero_cis',
	'snp_type', 'method_type', 'method')
dat <- rbind(dat_pwas[, cols_keep], dat_twas[, cols_keep])


#############################
## CLEAN UP DATA

# numeric pvalues
dat$PVAL <- as.numeric(dat$PVAL)

# FACTORS
dat$method <- factor(dat$method, 
	levels=c('cis-TWAS', 'cis/trans-TWAS', 
		'cis-PWAS', 'cis/trans-PWAS'))

# method type for later factoring
dat$method_type <- factor(dat$method_type, 
	levels=c('TWAS', 'PWAS'))

# snp type for later factoring
dat$snp_type <- factor(dat$snp_type, 
	levels=c('cis', 'cis/trans'))

# sort by CHROM, POS
dat <- dat[order(dat$CHROM, dat$GeneStart, dat$GeneEnd), ]
row.names(dat) <- NULL


#############################
# SIGNIFICANCE/SIGNIFICANT DATA

twas_sig_lvl <- 2.5e-6
pwas_sig_lvl <- 0.05 / 8922 # nrow(ct_pwas) = nrow(c_pwas) = 8922
## 5.604125e-06

# set significance level
dat$sig_level <- ifelse(dat$method_type == 'PWAS', pwas_sig_lvl, twas_sig_lvl)

# binary if significant
dat$sig <- as.integer(dat$adj_PVAL < dat$sig_level)

# significant dat
sig_dat <- dat[which(dat$sig == 1), ]


#############################
## SUMMARY TABLE

# summary table of results
res_dat <- cbind(
	# total number of genes
	tapply(dat, dat$method, nrow), 
	# number sig genes
	tapply(dat$sig, dat$method, function(x) {sum(x, na.rm=TRUE)}), 
	# number non-sig
	tapply(dat$sig, dat$method, function(x) {sum(as.integer(x == 0), na.rm=TRUE)}), 
	# number of NA genes
	tapply(dat$sig, dat$method, function(x) {sum(is.na(x))}), 
	# trans only genes
	tapply(as.integer(dat$zero_cis), dat$method, sum, na.rm=TRUE),
	# lambda_genomic_control
	tapply(dat$lambda, dat$method, unique)
	) |> as.data.frame()
colnames(res_dat) <- c('total_genes', 'sig_genes', 'nonsig_genes', 'NA_genes', 'trans_only')

# table for report
res_xtable <- xtable(res_dat,
	display = rep('d', 7)
	)
caption(res_xtable) <- 'xWAS Association Test Results Summary'
label(res_xtable) <- 'tab:summary'

print(res_xtable, 
	include.rownames = TRUE, 
	booktabs = TRUE, 
	table.placement = 'H', 
	sanitize.text.function = function(x){x},
	caption.placement = 'top',
	append = FALSE,
	file = here('tex', 'xtables.tex'))

cat('\\pagebreak\n\n', append = TRUE,
	file = here('tex', 'xtables.tex'))

#############################
## RISK GENES TABLES

sigtab <- sig_dat
sigtab$CHROM <- as.integer(sigtab$CHROM)
sigtab <- sigtab[, c('CHROM', 'GeneStart', 'GeneEnd', 'GeneName', 'method', 'adj_PVAL')]
colnames(sigtab) <- c('CHR', 'Start', 'End', 'Gene', 'method', 'p-value')
sigtab <- sigtab[order(sigtab$CHR, sigtab$`p-value`, sigtab$Start), ]

# make gene names italic
sigtab$Gene <-  paste0('\\textit{', sigtab$Gene, '}')

# print table for each method
for (method in unique(sigtab$method)){
	sigtab_method <- sigtab[sigtab$method == method, c('CHR', 'Start', 'End', 'Gene', 'p-value')]

	sigtab_xtable <- xtable(sigtab_method, 
		math.style.exponent=TRUE, display=c('d', 'd', 'd',  'd', 's', 'e'), include.rownames=FALSE)

	caption(sigtab_xtable) <- paste0('AD risk genes identified by ', method)
	label(sigtab_xtable) <- paste0('tab:', method, '_riskgenes')

	print(sigtab_xtable, 
		include.rownames = FALSE, 
		booktabs = TRUE, 
		table.placement = 'H', 
		tabular.environment = 'longtable', 
		sanitize.text.function = function(x){x},
		caption.placement = 'top',
		append = TRUE,
		file = here('tex', 'xtables.tex'))

	cat('\\pagebreak\n\n', append = TRUE,
		file = here('tex', 'xtables.tex'))
}

#############################
# VENN DIAGRAM

sig_dat$method <- factor(sig_dat$method, 
	levels=c('cis/trans-TWAS', 'cis-TWAS', 
		'cis-PWAS', 'cis/trans-PWAS'))
vennlist <- as.list(by(sig_dat, sig_dat$method, function(x) x[['GeneName']]))

p <- ggVennDiagram(vennlist, label_alpha=0, label='count') + 
	theme(legend.pos='None') + 
	scale_fill_gradient(low='white', high='#D92B26') +
	scale_x_continuous(expand=expansion(mult=0.11)) # expand so long labels fit
ggsave(here('figs', 'venn.pdf'), p)


#############################
# MANHATTAN PLOT

# dataframe for plots
mdat <- dat[complete.cases(dat),]

# add POS column
mdat$POS <- mdat$GeneStart

# add Pvalue column
mdat$Pvalue <- mdat$adj_PVAL

# factor in order
mdat$method <- factor(mdat$method, 
	levels=c('cis-TWAS', 'cis/trans-TWAS',
		'cis-PWAS', 'cis/trans-PWAS'))

# re-sort sig_data for upset_dat
sig_dat$method <- factor(sig_dat$method, 
	levels=c('cis-TWAS', 'cis/trans-TWAS', 'cis-PWAS', 'cis/trans-PWAS'))

upset_dat <- dcast(sig_dat, GeneName + CHROM + GeneStart + GeneEnd ~ method, sum, value.var='sig')

## get gene name lists for manplot
upset_dat$nsig <- rowSums(upset_dat[,c('cis/trans-TWAS', 'cis-TWAS', 'cis-PWAS', 'cis/trans-PWAS')])

# id'd by cis/trans-PWAS AND cis-PWAS
label_color_genes_cis_cistrans_pwas <- upset_dat[(upset_dat[['cis-PWAS']] == 1 & upset_dat[['cis/trans-PWAS']] == 1), 'GeneName']

# id'd by cis/trans-PWAS but not cis-PWAS
label_color_genes_cistrans_only_pwas <- upset_dat[(upset_dat[['cis-PWAS']] == 0 & upset_dat[['cis/trans-PWAS']] == 1), 'GeneName']

# id'd by cis/trans-TWAS AND cis-TWAS
label_color_genes_cis_cistrans_twas <- upset_dat[(upset_dat[['cis-TWAS']] == 1 & upset_dat[['cis/trans-TWAS']] == 1), 'GeneName']

# id'd by cis/trans-TWAS but not cis-TWAS
label_color_genes_cistrans_only_twas <- upset_dat[(upset_dat[['cis-TWAS']] == 0 & upset_dat[['cis/trans-TWAS']] == 1), 'GeneName']

id_genes <- unique(c(label_color_genes_cis_cistrans_pwas, label_color_genes_cistrans_only_pwas, label_color_genes_cis_cistrans_twas, label_color_genes_cistrans_only_twas))
mdat$label_text <- mdat$GeneName
mdat[with(mdat, !(GeneName %in% id_genes)), 'label_text'] <- ''

# don't label if sig for other method
mdat[with(mdat, (method_type=='PWAS') & !(GeneName %in% c(label_color_genes_cis_cistrans_pwas, label_color_genes_cistrans_only_pwas))), 'label_text'] <- ''
mdat[with(mdat, (method_type=='TWAS') & !(GeneName %in% c(label_color_genes_cis_cistrans_twas, label_color_genes_cistrans_only_twas))), 'label_text'] <- ''

## don't label for cis-cistrans overlap
mdat[with(mdat, (GeneName %in% label_color_genes_cis_cistrans_pwas) & (method_type == 'PWAS')), 'label_text'] <- ''
mdat[with(mdat, (GeneName %in% label_color_genes_cis_cistrans_twas) & (method_type == 'TWAS')), 'label_text'] <- ''

# label colors
mdat$label_color <- '#FD923F'

cis_cistrans_overlap_col <- '#D92B26'
mdat[with(mdat, (GeneName %in% label_color_genes_cis_cistrans_pwas) & (method_type == 'PWAS')), 'label_color'] <- cis_cistrans_overlap_col
mdat[with(mdat, (GeneName %in% label_color_genes_cis_cistrans_twas) & (method_type == 'TWAS')), 'label_color'] <- cis_cistrans_overlap_col

cistrans_only_col <- '#1032a9'
mdat[with(mdat, (GeneName %in% label_color_genes_cistrans_only_pwas) & (method_type == 'PWAS')), 'label_color'] <- cistrans_only_col
mdat[with(mdat, (GeneName %in% label_color_genes_cistrans_only_twas) & (method_type == 'TWAS')), 'label_color'] <- cistrans_only_col


# plot for slides
p_slides <- manhattan_plot(mdat, 
		max_iter=300000, 
		max_overlaps=100, 
		label_force=4, 
		theme=theme_gray(base_family = 'Fira Sans'), 
		plot_bg_col = '#FAFAFA', 
		strip_background=element_rect(color='#23373B', fill='#23373B')) +
		theme(
			strip.text = element_text(color='#FAFAFA', size=rel(0.9)), 
			strip.placement = 'outside') + 
	facet_grid(snp_type ~ method_type)
ggsave(here('figs', 'manplots_slides.png'), p_slides, height = 5, width = 8.5, dpi=400)


# plot for report
p_report <- manhattan_plot(mdat, 
		max_iter=300000, 
		max_overlaps=100, 
		label_force=4) + 
facet_grid(snp_type ~ method_type)
ggsave(here('figs', 'manplots_report.png'), p_report, height = 5, width = 8.5, dpi=400)

