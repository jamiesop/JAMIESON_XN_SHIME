

## Analysis & Figure scripts for SHIME project
## Purpose: xanthohumol's (XN) effect on gut microbiota metabolism in eubiotic and dysbiotic states
## Author: Paige Jamieson


## Environment ----------------------------


# Packages for Wrangling and Analysis
library(tidyverse)
library(DescTools) 
library(phyloseq)
library(vegan)

# Packages for figures
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggpubr)



## Load Data -----------------------------


# Meta Data
meta<- read.csv("~/Documents/SHIME project/SHIME/Data/SHIME_metadata.csv") %>% 
        mutate(across(reactor, ~ factor(.x, levels = c('AC', 'TC', 'DC'))))

# XN Metabolite Data
xn <- read.csv("~/Documents/SHIME project/SHIME/Data/SHIME_xn_metabs.csv") %>% 
        column_to_rownames("sample.id") %>% 
        rename_with(~paste0("XN_", .x)) %>% 
        rownames_to_column("sample.id")

# SCFA Data
scfa <- read.csv("~/Documents/SHIME project/SHIME/Data/SHIME_SCFA_Results_mM.csv") %>% 
        column_to_rownames("sample.id") %>% 
        rename_with(~paste0("FA_", .x)) %>% 
        rownames_to_column("sample.id")

# Join Sets
fulldata <- left_join(meta, xn) %>% 
        left_join(scfa) %>% 
        filter(period %in% c('control', 'treatment')) %>% 
        dplyr::filter(!day %in% c('30', '32'))


## Xanthohumol Metabolites -----------------------------------------

# Join data
xn_data <- fulldata %>% 
        pivot_longer(cols = c(starts_with('XN'), -XN_XN), names_to = 'metabolite', values_to = 'conc') %>% 
        mutate(across(metabolite, ~ factor(.x, 
                                           levels = c('XN_IXN', 'XN_DXN', 'XN_x8PN', 'XN_DDXN', 'XN_x6PN')))) %>% 
        dplyr::filter(week != 3) %>% 
        mutate_at('week', ~ ifelse(week == '4', 'B', 
                                   ifelse(week == '5', 'T1',
                                          ifelse(week == '6', 'T2', 'T3')))) %>% 
        mutate_at('reactor', ~ factor(.x, levels = c('AC', 'TC', 'DC')))

# Summary stats
xn_summ <- xn_data %>% 
        group_by(SHIME, reactor, week, metabolite) %>% 
        summarize(mConc = mean(conc),
                  serConc = ser(conc))

# color palette for XN Metabolites
xn_pal <- c( '#B6D886', '#FD8B8B', '#F56147', '#F7C059',  '#529163')
names(xn_pal) <- c('XN_IXN', 'XN_DXN', 'XN_x8PN', 'XN_DDXN', 'XN_x6PN')

# H-SHIME
shimeA <- xn_summ %>% 
        dplyr::filter(SHIME == 'A') %>% 
        dplyr::filter(week != 'B')

# D-SHIME
shimeB <- xn_summ %>% 
        dplyr::filter(SHIME == 'B') %>% 
        dplyr::filter(week != 'B')


# Stacked Bar Charts -------------------------

# Stacked bar chart of XN Metabolties (SHIME A)
barA <- ggplot(shimeA, aes(week, mConc, fill = metabolite, group = SHIME)) +
        geom_bar(stat = 'identity') + 
        facet_wrap(~ reactor, ncol = 1) +
        scale_fill_manual(values = xn_pal, 
                          labels = c('IXN', 'DXN', '8PN', 'DDXN', '6PN'), 
                          name = 'Metabolite',
                          guide = 'none') +
        theme_bw() +
        scale_y_continuous(limits = c(0, 135)) + 
        labs(x = '', y = '') +
        theme(axis.title=element_text(size=24),
              axis.text.x=element_text(size=22),
              axis.text.y = element_text(size=0),
              axis.ticks.y = element_blank(),
              strip.text = element_text(size=24)) +
        theme(strip.text.y = element_text(size = 20, color = 'black'),
              strip.text.x = element_text(size = 20, color = 'black'),
              strip.background = element_rect(color =  'white', fill = 'white'),
        )

# Stacked bar chart of XN Metabolties (SHIME B)
barB <- ggplot(shimeB, aes(week, mConc, fill = metabolite, group = SHIME)) +
        geom_bar(stat = 'identity') + 
        facet_wrap(~ reactor, ncol = 1) +
        scale_fill_manual(values = xn_pal, 
                          labels = c('IXN', 'DXN', '8PN', 'DDXN', '6PN'), 
                          name = 'Metabolite',
                          guide = 'none') +
        theme_bw() +
        scale_y_continuous(limits = c(0, 135)) + 
        labs(x = '', y = '') +
        theme(axis.title=element_text(size=24),
              axis.text=element_text(size=22),
              axis.text.y = element_text(size=0),
              axis.ticks.y = element_blank(),
              strip.text = element_text(size=24)) +
        theme(strip.text.y = element_text(size = 20, color = 'black'),
              strip.text.x = element_text(size = 20, color = 'black'),
              strip.background = element_rect(color =  'white', fill = 'white'),
        )


# Longitudinal Plots ----------------------

# H-SHIME
xn_shimeA <- xn_summ %>% 
        dplyr::filter(SHIME == 'A')
# D-SHIME
xn_shimeB <- xn_summ %>% 
        dplyr::filter(SHIME == 'B')

# shime A [H-SHIME]
longA <- ggplot(xn_shimeA, aes(week, mConc, fill = metabolite, group = metabolite, color = metabolite)) + 
        geom_point(size = 2, stroke = 1.5, show.legend = F) +
        geom_path(size = 1.5, show.legend = F) +
        geom_errorbar(aes(ymin = mConc - serConc, ymax = mConc + serConc), width = 0.2, size = 1.5) +
        facet_wrap(~reactor, ncol = 1) + 
        scale_color_manual(values = xn_pal, 
                           labels = c('IXN', 'DXN', '8PN', 'DDXN', '6PN'), 
                           name = 'Metabolite',
                           guide = 'none') +
        scale_fill_manual(values = xn_pal) +
        theme_bw() +
        scale_y_continuous(limits = c(0, 135)) +
        labs(x = '', y = 'nM') +
        theme(axis.title=element_text(size=24),
              axis.text=element_text(size=22),
              strip.text = element_text(size=24)) +
        theme(strip.text.y = element_text(size = 20, color = 'black'),
              strip.text.x = element_text(size = 20, color = 'black'),
              strip.background = element_rect(color =  'white', fill = 'white'),
        )


# shime B [D-SHIME]
longB <- ggplot(xn_shimeB, aes(week, mConc, fill = metabolite, group = metabolite, color = metabolite)) + 
        geom_point(size = 2, stroke = 1.5, show.legend = F) +
        geom_path(size = 1.5, show.legend = F) +
        geom_errorbar(aes(ymin = mConc - serConc, ymax = mConc + serConc), width = 0.2, size = 1.5) +
        facet_wrap(~reactor, ncol = 1) + 
        scale_color_manual(values = xn_pal, 
                           labels = c('IXN', 'DXN', '8PN', 'DDXN', '6PN'), 
                           name = 'Metabolite',
                           guide = 'none') +
        scale_fill_manual(values = xn_pal) +
        theme_bw() +
        scale_y_continuous(limits = c(0, 135)) +
        labs(x = '', y = 'nM') +
        theme(axis.title=element_text(size=24),
              axis.text=element_text(size=22),
              strip.text = element_text(size=24)) +
        theme(strip.text.y = element_text(size = 20, color = 'black'),
              strip.text.x = element_text(size = 20, color = 'black'),
              strip.background = element_rect(color =  'white', fill = 'white'),
        )

# Extract legend
L <- ggplot(xn_shimeA, aes(week, mConc, fill = metabolite, group = SHIME)) +
        geom_bar(stat = 'identity') + 
        facet_wrap(~ reactor, ncol = 1) +
        scale_fill_manual(values = xn_pal, 
                          labels = c('IXN', 'DXN', '8PN', 'DDXN', '6PN'), 
                          name = 'Metabolite') +
        theme_bw() +
        scale_y_continuous(limits = c(0, 125)) + 
        labs(x = 'week', y = '') +
        theme(axis.title=element_text(size=14),
              axis.text=element_text(size=12),
              strip.text = element_text(size=14)) +
        theme(legend.key.size = unit(2, 'cm'), #change legend key size
              legend.key.height = unit(2, 'cm'), #change legend key height
              legend.key.width = unit(2, 'cm'), #change legend key width
              legend.title = element_text(size=30), #change legend title font size
              legend.text = element_text(size=26)) #change legend text font size

leg <- get_legend(L)

# [SAVE 700w x 600h]
leg <- as_ggplot(leg) 

# XN Metabolites Longitudinal and Stacked Bar Plots
# [SAVE 1200w x 800h]
plot_grid(longA, barA, longB, barB, nrow = 1, rel_widths = c(2,1,2,1))


## Short Chain Fatty Acids -----------------------------------------

# One-Way ANOVA
scfa_aov <- fulldata %>% 
        dplyr::select(SHIME, reactor, week, FA_Acetic, FA_Propionic, FA_Butyric) %>% 
        pivot_longer(cols = c("FA_Acetic", 'FA_Propionic', 'FA_Butyric'), 
                     names_to = 'measure', values_to = 'value') %>% 
        group_by(SHIME, reactor, measure) %>% 
        nest() %>% 
        mutate(mod = purrr::map(data, function(x) aov(value ~ week, data = x))) %>% 
        mutate(pval = purrr::map_dbl(mod, function(x) summary(x)[[1]][[5]][[1]] )) %>% 
        mutate(incl = purrr::map_chr(pval, function(x) ifelse(x <= 0.05, 'Y', 'N')))


# Post Hoc Dunnett's Test with Bonferroni correction
scfa_dun <- scfa_aov %>% 
        mutate(dunmod = purrr::map(data, function(x) DunnettTest(x$value, x$week, data = x))) %>% 
        mutate(pval_t1 = purrr::map_dbl(dunmod, function(x) x$`3`[[2,4]])) %>%
        mutate(pval_t2 = purrr::map_dbl(dunmod, function(x) x$`3`[[3,4]])) %>%
        mutate(pval_t3 = purrr::map_dbl(dunmod, function(x) x$`3`[[4,4]])) %>% 
        ungroup() %>% 
        filter(incl == 'Y') %>% 
        mutate(across(starts_with('pval_'), ~ p.adjust(.x, method = 'bonferroni')))

# Reorder for figure
scfa4plot <- fulldata %>% 
        mutate(time = ifelse(week == 3, 'B', 
                             ifelse(week == 5, 'T1', 
                                    ifelse(week == 6, 'T2', 'T3')))) %>% 
        dplyr::select(SHIME, reactor, time, FA_Acetic, FA_Propionic, FA_Butyric) %>% 
        pivot_longer(cols = c("FA_Acetic", 'FA_Propionic', 'FA_Butyric'), 
                     names_to = 'measure', values_to = 'value') %>% 
        group_by(SHIME, reactor, time, measure) %>%
        summarise(mFA = mean(value),
                  sdFA = sd(value)) %>% 
        mutate_at('measure', ~ factor(.x, levels = c('FA_Acetic', 'FA_Propionic', 'FA_Butyric')))

# Significance annotations
scfa_sigs <- scfa_dun %>% 
        dplyr::select(-data, -mod, -pval) %>% 
        pivot_longer(starts_with('pval'), names_to = 'time', values_to = 'pval') %>% 
        mutate(time = ifelse(time == 'pval_t1', 'T1', ifelse(time == 'pval_t2', 'T2', 'T3'))) %>% 
        filter(pval <= 0.05) %>% 
        filter(incl == 'Y') %>% 
        left_join(scfa4plot) %>% 
        mutate(mFA_new = mFA + sdFA + 0.9) %>% 
        dplyr::select(-mFA, -sdFA, - incl) %>% 
        rename(mFA = mFA_new) %>% 
        mutate_at('measure', ~ factor(.x, levels = c('FA_Acetic', 'FA_Propionic', 'FA_Butyric'))) %>% 
        mutate_at('SHIME', ~ factor(.x, levels = c('A', 'B')))

# Significance annotation
ann_text1 <- scfa_sigs %>% 
        filter(pval <= 0.05 & pval > 0.01)

ann_text2A <- scfa_sigs %>% 
        filter(SHIME == 'A') %>% 
        filter(pval <= 0.01)    

ann_text2B <- scfa_sigs %>% 
        filter(SHIME == 'B') %>% 
        filter(pval <= 0.01)  

scfalabs <- c('Acetic acid', 'Propionic acid', 'Butyric acid')
names(scfalabs) <- c('FA_Acetic', 'FA_Propionic', 'FA_Butyric')

# SCFA plot
# [Save 1200w x 1000h]
ggplot(scfa4plot, aes(time, mFA, group = SHIME, shape = SHIME)) + 
        geom_point(size = 2, stroke = 1.5, position = position_dodge(width = 0.2)) +
        geom_path(size = 1, position = position_dodge(width = 0.2)) +
        geom_errorbar(aes(ymin = mFA - sdFA, ymax = mFA + sdFA), width = 0.15, size = 1, position = position_dodge(width = 0.2)) +
        scale_shape_manual(values = c(0, 2), labels = c('H-SHIME', 'D-SHIME')) +
        facet_grid(measure ~ reactor, labeller = labeller(measure = scfalabs)) +
        theme_bw() +
        theme(panel.border = element_rect(size = 1),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 20),
              axis.ticks = element_blank(),
              legend.text = element_text(size = 18)) +
        theme(strip.text.y = element_text(size = 18, color = 'black'),
              strip.text.x = element_text(size = 20, color = 'black'),
              strip.background = element_rect(color =  'white', fill = 'white')) +
        geom_text(data = ann_text1, aes(time, mFA, group = SHIME), nudge_x = -0.05, label = "*", size = 9, color = 'black') +
        geom_text(data = ann_text2B, aes(time, mFA, group = SHIME), nudge_x = 0.05, label = "**", size = 9, color = 'black') +
        geom_text(data = ann_text2A, aes(time, mFA, group = SHIME), nudge_x = -0.05, label = "**", size = 9, color = 'black') +
        labs(x = "", y = 'mM', shape = "")



## Colonic Microbiome -----------------------------------------------

# Load phyloseq object
shime_ps <- readRDS("~/Documents/16S/SHIME/shime_ps.rds")


# Create placeholder names for un-annotated genera using family
renames <- rownames(tax_table(shime_ps)[is.na(tax_table(shime_ps)[, 'Genus'])])
taxdf <- tax_table(shime_ps)[renames,]
renamed_genus <- unname(sapply(taxa_names(taxdf), function(x) paste0('f_', taxdf[x, 'Family'])))
tax_table(shime_ps)[renames, 'Genus'] <- renamed_genus

rename_order <- tax_table(shime_ps) %>% data.frame() %>% filter(., grepl("f_NA", Genus)) %>% rownames()
taxadf <- tax_table(shime_ps)[rename_order,]
renamed_from_order <- unname(sapply(taxa_names(taxadf), function(x) paste0('o_', taxadf[x, 'Order'], '_')))
tax_table(shime_ps)[rename_order, 'Genus'] <- renamed_from_order


# Raw count data - Remove stabilization and wash-out period
shime_ps_filt <- shime_ps %>% subset_samples(period %in% c('control', 'treatment'))


# Agglomerate to Genus level
ps_genera <- shime_ps_filt %>% tax_glom(taxrank = 'Genus') #ASVs = 187

# Filter out taxa seen fewer than 3 times in less than 20% of samples
gen_counts <- ps_genera %>% filter_taxa(function(x) sum(x > 3) > (0.2*length(x)), TRUE) # ASVs: 88

#Filter out low abundance (>1e-5) taxa >>> does not change number of ASVs 
gen_filtered <- gen_counts %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE) # ASVs = 88

# Transform to relative abundance
gen_relab <- gen_filtered %>% transform_sample_counts(function(x) x / sum(x))


# Taxa table
tax <- tax_table(shime_ps) %>% 
        data.frame() %>% 
        rownames_to_column('ASV')


## Alpha Diversity --------------------------------------------------


# Meta data        
meta <- shime_ps_filt %>% sample_data() %>% data.frame() %>% rownames_to_column('sample.id')

# Alpha diversity measures
adiv <- shime_ps_filt %>%
        estimate_richness(measures = c('Observed', 'Shannon', 'Simpson')) %>% 
        rownames_to_column('sample.id') %>% 
        left_join(., meta) %>% 
        filter(week != 3) %>% 
        dplyr::filter(!day %in% c('30', '32')) %>% 
        mutate_at('reactor', ~ factor(.x, levels = c('AC', 'TC', 'DC'))) %>% 
        mutate_at('SHIME', ~ factor(.x, levels = c('A', 'B')))

# One-Way ANOVA
adiv_aov <- adiv %>% 
        dplyr::select(SHIME, reactor, week, Observed, Shannon) %>% 
        pivot_longer(cols = c('Observed', 'Shannon'), names_to = 'measure', values_to = 'value') %>% 
        group_by(SHIME, reactor, measure) %>%
        nest() %>% 
        mutate(mod = purrr::map(data, function(x) aov(value ~ week, data = x))) %>% 
        mutate(pval = purrr::map_dbl(mod, function(x) summary(x)[[1]][[5]][[1]] )) %>% 
        mutate(incl = purrr::map_chr(pval, function(x) ifelse(x <= 0.05, 'Y', 'N')))


# Post Hoc Dunnett's Test with Bonferroni correction
adiv_dun <- adiv_aov %>% 
        mutate(dunmod = purrr::map(data, function(x) DunnettTest(x$value, x$week, data = x))) %>% 
        mutate(pval_t1 = purrr::map_dbl(dunmod, function(x) x$`4`[[1,4]])) %>%
        mutate(pval_t2 = purrr::map_dbl(dunmod, function(x) x$`4`[[2,4]])) %>%
        mutate(pval_t3 = purrr::map_dbl(dunmod, function(x) x$`4`[[3,4]])) %>% 
        ungroup() %>%
        dplyr::filter(incl == 'Y') %>% 
        mutate(across(starts_with('pval_'), ~ p.adjust(.x, method = 'bonferroni')))


# Reorder for figure
adiv4plot <- adiv %>% 
        filter(week != 3) %>% 
        mutate(time = ifelse(week == 4, 'B', 
                             ifelse(week == 5, 'T1', 
                                    ifelse(week == 6, 'T2', 'T3')))) %>% 
        dplyr::select(SHIME, reactor, time, Observed, Shannon) %>% 
        pivot_longer(cols = c("Observed", 'Shannon'), 
                     names_to = 'measure', values_to = 'value') %>% 
        group_by(SHIME, reactor, time, measure) %>%
        summarise(mFA = mean(value),
                  sdFA = sd(value)) %>% 
        mutate_at('measure', ~ factor(.x, levels = c('Observed', 'Shannon'))) %>% 
        mutate_at('reactor', ~ factor(.x, levels = c('AC', 'TC', 'DC')))


# For sig annotations
adiv_stat_ann <- adiv_dun %>% 
        dplyr::select(-mod, -pval, -dunmod, -data) %>% 
        pivot_longer(starts_with('pval'), names_to = 'time', values_to = 'pval') %>% 
        mutate(time = ifelse(time == 'pval_t1', 'T1', ifelse(time == 'pval_t2', 'T2', 'T3'))) %>% 
        dplyr::filter(pval <= 0.05) %>% 
        left_join(adiv4plot) %>% 
        mutate(mFA_new = mFA + sdFA + 10) %>% 
        dplyr::select(-mFA, -sdFA) %>% 
        rename(mFA = mFA_new) 

# Significance annotation
ann_below.05.Obs <- adiv_stat_ann %>% 
        filter(measure == 'Observed') %>% 
        filter(pval <= 0.05 & pval > 0.01)

ann_below.05.Shan <- adiv_stat_ann %>% 
        filter(measure == 'Shannon') %>% 
        mutate_at('mFA', ~ .x - 9.9) %>% 
        filter(pval <= 0.05 & pval > 0.01)
        
ann_below.01 <- adiv_stat_ann %>% 
        filter(pval <= 0.01)

# Alpha Diversity plot
# [SAVE 1200w x 650h]
ggplot(adiv4plot, aes(time, mFA, group = SHIME, shape = SHIME)) +
        geom_point(size = 2, stroke = 1.5, position = position_dodge(width = 0.2)) +
        geom_path(size = 1, position = position_dodge(width = 0.2)) +
        geom_errorbar(aes(ymin = mFA - sdFA, ymax = mFA + sdFA), width = 0.15, size = 1, position = position_dodge(width = 0.2)) +
        facet_grid(measure ~ reactor, scales = 'free') +
        scale_shape_manual(values = c(0, 2), labels = c('H-SHIME', 'D-SHIME')) +
        theme_bw() +
        theme(panel.border = element_rect(size = 1),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 20),
              axis.ticks = element_blank(),
              legend.text = element_text(size = 18)) +
        theme(strip.text.y = element_text(size = 18, color = 'black'),
              strip.text.x = element_text(size = 20, color = 'black'),
              strip.background = element_rect(color =  'white', fill = 'white')) +
        labs(x = "", y = '', shape = '') +
        geom_text(data = ann_below.05.Shan, aes (time, mFA, group = SHIME), nudge_x = -0.05, label = '*', size = 9, color = 'black') +
        geom_text(data = ann_below.05.Obs, aes (time, mFA, group = SHIME), nudge_x = 0.05, label = '*', size = 9, color = 'black') +
        geom_text(data = ann_below.01, aes (time, mFA, group = SHIME), nudge_x = -0.05, label = '**', size = 9, color = 'black')



# Differential Abundance ------------------------------------------------

# CLR Transform counts
gen_clr <- gen_filtered %>% 
        # rarefy to even depth
        rarefy_even_depth(rngseed = 44) %>% 
        otu_table() %>% 
        data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1) %>% 
        rownames_to_column('sample.id') # still 88 ASVs

#Pull metadata 
metadata <- gen_filtered %>%  
        sample_data() %>% 
        data.frame() %>%
        rownames_to_column('sample.id') 

# Join data
micro <- left_join(metadata, gen_clr) %>% 
        mutate(across(week, ~ factor(.x)))

# One-way ANOVA
clr_aov <- micro %>% 
        filter(day != '29' & week != '4') %>% 
        dplyr::select(SHIME, reactor, period, week, sampreact, starts_with('ASV')) %>% 
        pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = 'clr') %>% 
        group_by(SHIME, reactor, ASV) %>% 
        nest() %>% 
        mutate(mod = purrr::map(data, function(x) aov(clr ~ week, data = x))) %>% 
        mutate(pval = purrr::map_dbl(mod, function(x) summary(x)[[1]][[5]][[1]])) %>% 
        mutate(incl = purrr::map_chr(pval, function(x) ifelse( x <= 0.05, 'Y', 'N')))

# Post Hoc Dunnett's Test with bonferroni
clr_dun <- clr_aov %>% 
        mutate(dunmod = purrr::map(data, function(x) DunnettTest(x$clr, x$week, data = x))) %>% 
        mutate(pval_t1 = purrr::map_dbl(dunmod, function(x) x$`3`[[1,4]])) %>%
        mutate(pval_t2 = purrr::map_dbl(dunmod, function(x) x$`3`[[2,4]])) %>%
        mutate(pval_t3 = purrr::map_dbl(dunmod, function(x) x$`3`[[3,4]])) %>% 
        ungroup() %>%
        dplyr::filter(incl == 'Y') %>% 
        mutate(across(starts_with('pval_'), ~ p.adjust(.x, method = 'bonferroni'))) %>% 
        left_join(tax) %>% 
        mutate(genspec = ifelse(is.na(Species), Genus, paste0(Genus, "\n", Species)))

# Relative Abundance Table
relab <- gen_relab %>% 
        otu_table() %>% 
        data.frame() %>% 
        rownames_to_column('sample.id')

# Significant features
asv_sig_names <- clr_dun %>% 
        dplyr::filter(pval_t1 <= 0.05 & pval_t1 > 0 | 
                              pval_t2 <= 0.05 & pval_t2 > 0 | 
                              pval_t3 <= 0.05 & pval_t3 > 0) %>% 
        pull(ASV)

# Relative abundance data reordered for figure
micro_relab <- gen_relab %>% 
        sample_data() %>% 
        data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        left_join(relab) %>% 
        dplyr::select(sample.id, SHIME, reactor, period, day, week, all_of(asv_sig_names)) %>% 
        filter(day != '29' & week != '4') %>% 
        pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = 'relab') %>% 
        mutate(week = ifelse(week == 3, 'B', 
                             ifelse(week == 5, 'T1', 
                                    ifelse(week == 6, 'T2', 'T3')))) %>% 
        mutate(across('reactor', ~factor(.x, levels = c('AC', 'TC', 'DC'), ordered = T)))

# Summary stats
micro_summ <- micro_relab %>% 
        group_by(SHIME, reactor, week, ASV) %>% 
        summarize(mRelab = mean(relab),
                  sdRelab = sd(relab))

# Significant featuers for H-SHIME
asv_sig_namesA <- clr_dun %>% 
        dplyr::filter(pval_t1 <= 0.05 & pval_t1 > 0 | 
                              pval_t2 <= 0.05 & pval_t2 > 0 | 
                              pval_t3 <= 0.05 & pval_t3 > 0) %>% 
        filter(SHIME == 'A') 

# Significant features for D-SHIME
asv_sig_namesB <- clr_dun %>% 
        dplyr::filter(pval_t1 <= 0.05 & pval_t1 > 0 | 
                              pval_t2 <= 0.05 & pval_t2 > 0 | 
                              pval_t3 <= 0.05 & pval_t3 > 0) %>% 
        filter(SHIME == 'B') 

namesA <- asv_sig_namesA$genspec
names(namesA) <- asv_sig_namesA$ASV

namesB <- asv_sig_namesB$genspec
names(namesB) <- asv_sig_namesB$ASV


# Sig annotations (H-SHIME)
annoA <- asv_sig_namesA %>% 
        dplyr::select(SHIME, reactor, ASV, genspec, pval_t1, pval_t2, pval_t3) %>% 
        pivot_longer(cols = starts_with('pval_'), names_to = 'week', values_to = 'value') %>% 
        mutate(week = ifelse(week == 'pval_t1', 'T1', ifelse(week == 'pval_t2', 'T2', 'T3'))) %>% 
        dplyr::filter(value <= 0.05) %>% 
        left_join(micro_summ %>% dplyr::select(SHIME, reactor, week, ASV, mRelab, sdRelab)) %>% 
        mutate(mRelab = mRelab + sdRelab + 0.002)

# Sig annotations < 0.05 (D-SHIME)
annoB.below.05 <- asv_sig_namesB %>% 
        dplyr::select(SHIME, reactor, ASV, genspec, pval_t1, pval_t2, pval_t3) %>% 
        pivot_longer(cols = starts_with('pval_'), names_to = 'week', values_to = 'value') %>% 
        mutate(week = ifelse(week == 'pval_t1', 'T1', ifelse(week == 'pval_t2', 'T2', 'T3'))) %>% 
        dplyr::filter(value <= 0.05 & value > 0.01) %>% 
        left_join(micro_summ %>% dplyr::select(SHIME, reactor, week, ASV, mRelab, sdRelab)) %>% 
        mutate(mRelab = mRelab + sdRelab + 0.001) %>% 
        dplyr::filter(!is.na(mRelab))

# Sig annotations < 0.01 (D-SHIME)
annoB.below.01 <- asv_sig_namesB %>% 
        dplyr::select(SHIME, reactor, ASV, genspec, pval_t1, pval_t2, pval_t3) %>% 
        pivot_longer(cols = starts_with('pval_'), names_to = 'week', values_to = 'value') %>% 
        mutate(week = ifelse(week == 'pval_t1', 'T1', ifelse(week == 'pval_t2', 'T2', 'T3'))) %>% 
        dplyr::filter(value <= 0.01 & value > 0.001) %>% 
        left_join(micro_summ %>% dplyr::select(SHIME, reactor, week, ASV, mRelab, sdRelab)) %>% 
        mutate(mRelab = mRelab + sdRelab + 0.0002) %>% 
        dplyr::filter(!is.na(mRelab))

# Sig annotations < 0.001 (D-SHIME)
annoB.below.001.AG <- asv_sig_namesB %>% 
        filter(ASV == 'ASV13') %>% 
        dplyr::select(SHIME, reactor, ASV, genspec, pval_t1, pval_t2, pval_t3) %>% 
        pivot_longer(cols = starts_with('pval_'), names_to = 'week', values_to = 'value') %>% 
        mutate(week = ifelse(week == 'pval_t1', 'T1', ifelse(week == 'pval_t2', 'T2', 'T3'))) %>% 
        dplyr::filter(value <= 0.001 & value > 0) %>% 
        left_join(micro_summ %>% dplyr::select(SHIME, reactor, week, ASV, mRelab, sdRelab)) %>% 
        mutate(mRelab = mRelab + sdRelab + 0.0006) %>% 
        dplyr::filter(!is.na(mRelab))

# Sig annotations < 0.001 (D-SHIME)    
annoB.below.001.NS <- asv_sig_namesB %>% 
        filter(ASV == 'ASV167') %>% 
        dplyr::select(SHIME, reactor, ASV, genspec, pval_t1, pval_t2, pval_t3) %>% 
        pivot_longer(cols = starts_with('pval_'), names_to = 'week', values_to = 'value') %>% 
        mutate(week = ifelse(week == 'pval_t1', 'T1', ifelse(week == 'pval_t2', 'T2', 'T3'))) %>% 
        dplyr::filter(value <= 0.001 & value > 0) %>% 
        left_join(micro_summ %>% dplyr::select(SHIME, reactor, week, ASV, mRelab, sdRelab)) %>% 
        mutate(mRelab = mRelab + sdRelab + 0.00002) %>% 
        dplyr::filter(!is.na(mRelab))

# sig ASV names    
microA <- micro_summ %>% 
        dplyr::filter(SHIME == 'A') %>% 
        dplyr::filter(ASV %in% annoA$ASV)

# sig ASV names
microB <- micro_summ %>% 
        dplyr::filter(SHIME == 'B') %>% 
        dplyr::filter(ASV %in% c(annoB.below.001$ASV, annoB.below.01$ASV, annoB.below.05$ASV))



# Plot differentiall abundance taxa in H-SHIME 
plotA <- ggplot(microA, aes(week, mRelab, group = SHIME, shape = SHIME)) +
        geom_point(size = 2, stroke = 1.5, position = position_dodge(width = 0.2), show.legend = F) +
        geom_path(size = 1, position = position_dodge(width = 0.2)) +
        geom_errorbar(aes(ymin = mRelab - sdRelab, ymax = mRelab + sdRelab), width = 0.15, size = 1, position = position_dodge(width = 0.2)) +
        facet_grid(~ASV ~ factor(reactor), scales = 'free_y', labeller = labeller(ASV = namesA)) +
        #scale_shape_manual(values = c(0, 2), labels = c('H-SHIME', 'D-SHIME')) +
        theme_bw() +
        theme(panel.border = element_rect(size = 1),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 20),
              axis.ticks = element_blank(),
              legend.text = element_text(size = 18)) +
        theme(strip.text.y = element_text(size = 18, color = 'black'),
              strip.text.x = element_text(size = 20, color = 'black'),
              strip.background = element_rect(color =  'white', fill = 'white')) +
        labs(x = "", y = '', shape = '') +
        geom_text(data = annoA, 
                  aes (week, mRelab, group = SHIME), 
                  label = '*', 
                  size = 9, color = 'black')
# SHIME labels
shimelabs <- c('H-SHIME', 'D-SHIME')
names(shimelabs) <- c('A', 'B')

# Plot to extract SHIME title
plotAtitle <- ggplot(microA, aes(week, mRelab, group = SHIME, shape = SHIME)) +
        geom_point(size = 2, stroke = 1.5, position = position_dodge(width = 0.2), show.legend = F) +
        geom_path(size = 1, position = position_dodge(width = 0.2)) +
        geom_errorbar(aes(ymin = mRelab - sdRelab, ymax = mRelab + sdRelab), width = 0.15, size = 1, position = position_dodge(width = 0.2)) +
        facet_grid(~SHIME ~ factor(reactor), scales = 'free_y', labeller = labeller(SHIME = shimelabs)) +
        #scale_shape_manual(values = c(0, 2), labels = c('H-SHIME', 'D-SHIME')) +
        theme_bw() +
        theme(panel.border = element_rect(size = 1),
              axis.title = element_text(size = 24),
              axis.text = element_text(size = 20),
              axis.ticks = element_blank(),
              legend.text = element_text(size = 18)) +
        theme(strip.text.y = element_text(size = 28, color = 'black'),
              strip.text.x = element_text(size = 20, color = 'black'),
              strip.background = element_rect(color =  'white', fill = 'white')) +
        labs(x = "", y = 'Relative Abundance', shape = '') +
        geom_text(data = annoA, 
                  aes (week, mRelab, group = SHIME), 
                  label = '**', 
                  size = 9, color = 'black')

# Extract SHIME title
Atitle <- cowplot::get_plot_component(plotAtitle, "strip-r-1", return_all = T)
fullplotA <- cowplot::plot_grid(NULL, plotA, NULL, Atitle, ncol = 4, rel_widths = c(1, 20, 0.5, 2))

relablab <- cowplot::get_plot_component(plotAtitle, "ylab-l", return_all = T)
relablab <- cowplot::ggdraw(relablab)
 
# Plot differentiall abundance taxa in D-SHIME        
plotB <- ggplot(microB, aes(week, mRelab, group = SHIME, shape = SHIME)) +
        geom_point(size = 2, stroke = 1.5, position = position_dodge(width = 0.2), show.legend = F) +
        geom_path(size = 1, position = position_dodge(width = 0.2)) +
        geom_errorbar(aes(ymin = mRelab - sdRelab, ymax = mRelab + sdRelab), width = 0.15, size = 1, position = position_dodge(width = 0.2)) +
        #facet_grid(~Genus+Species ~ factor(reactor), scales = 'free', labeller = label_wrap_gen(multi_line = TRUE)) +
        facet_grid(~ASV ~ factor(reactor), scales = 'free', labeller = labeller(ASV = namesB)) +
        #scale_shape_manual(values = c(0, 2), labels = c('H-SHIME', 'D-SHIME')) +
        theme_bw() +
        theme(panel.border = element_rect(size = 1),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 20),
              axis.ticks = element_blank(),
              legend.text = element_text(size = 18)) +
        theme(strip.text.y = element_text(size = 18, color = 'black'),
              strip.text.x = element_text(size = 20, color = 'black'),
              strip.background = element_rect(color =  'white', fill = 'white')) +
        labs(x = "", y = '', shape = '') +
        geom_text(data = annoB.below.05, 
                  aes (week, mRelab, group = SHIME), 
                  label = '*', 
                  size = 9, color = 'black') +
        geom_text(data = annoB.below.01, 
                  aes (week, mRelab, group = SHIME), 
                  label = '**', 
                  size = 9, color = 'black') +
        geom_text(data = annoB.below.001.AG, 
                  aes (week, mRelab, group = SHIME), 
                  label = '***', 
                  size = 9, color = 'black') +
        geom_text(data = annoB.below.001.NS, 
                  aes (week, mRelab, group = SHIME), 
                  label = '***', 
                  size = 9, color = 'black')

# Plot to extract SHIME title
plotBtitle <- ggplot(microB, aes(week, mRelab, group = SHIME, shape = SHIME)) +
        geom_point(size = 2, stroke = 1.5, position = position_dodge(width = 0.2), show.legend = F) +
        geom_path(size = 1, position = position_dodge(width = 0.2)) +
        geom_errorbar(aes(ymin = mRelab - sdRelab, ymax = mRelab + sdRelab), width = 0.15, size = 1, position = position_dodge(width = 0.2)) +
        facet_grid(~SHIME ~ factor(reactor), scales = 'free', labeller = labeller(SHIME = shimelabs)) +
        #scale_shape_manual(values = c(0, 2), labels = c('H-SHIME', 'D-SHIME')) +
        theme_bw() +
        theme(panel.border = element_rect(size = 1),
              axis.title = element_text(size = 20),
              axis.text = element_text(size = 20),
              axis.ticks = element_blank(),
              legend.text = element_text(size = 18)) +
        theme(strip.text.y = element_text(size = 28, color = 'black'),
              strip.text.x = element_text(size = 20, color = 'black'),
              strip.background = element_rect(color =  'white', fill = 'white')) +
        labs(x = "", y = '', shape = '') +
        geom_text(data = annoB.below.05, 
                  aes (week, mRelab, group = SHIME), 
                  label = '*', 
                  size = 9, color = 'black') +
        geom_text(data = annoB.below.01, 
                  aes (week, mRelab, group = SHIME), 
                  label = '**', 
                  size = 9, color = 'black') +
        geom_text(data = annoB.below.001, 
                  aes (week, mRelab, group = SHIME), 
                  label = '***', 
                  size = 9, color = 'black')

# Extract SHIME title
Btitle <- cowplot::get_plot_component(plotBtitle, "strip-r-1", return_all = T)
fullplotB <- plot_grid(plotB, Btitle, ncol = 2, rel_widths = c(10,1))

# [SAVE 1200w x 1700h]
fullplot1 <- plot_grid(fullplotA, fullplotB, ncol = 1, rel_heights = c(1.6,7))
fullplot <- plot_grid(relablab, fullplot1, ncol = 2, rel_widths = c(1,30))



# METABOLOMICS -----------------------------------------------------------------


# Read in Progenesis QI datasets

# Positive ion mode
pos <- read.csv("~/Documents/SHIME Project/SHIME/Metabolomics/SHIME_POS_092424.csv", 
                skip = 2, na.string = c(""))

# Negative ion mode
neg <- read.csv("~/Documents/SHIME Project/SHIME/Metabolomics/SHIME_NEG_092424.csv", 
                skip = 2, na.strings = c(""))

# Read in sample identifiers
posids <- read.csv("~/Documents/SHIME Project/SHIME/Data/pos_sample_identifiers.csv")
negids <- posids %>% mutate(run.id = gsub("pos", 'neg', run.id))

# Read in meta data
meta <- read.csv("~/Documents/SHIME Project/SHIME/Data/SHIME_metadata.csv") %>% 
        mutate_at('reactor', ~factor(.x, levels = c('AC', 'TC', 'DC')),
                  'period', ~factor(.x, levels = c('stabilization', 'control', 'treatment', 'wash out')))

# Positive ion mode Compound IDs
posIDs <- pos %>%
        dplyr::select(!matches(c("^X07032024.*\\d\\.1$", "^X07032024.*\\d"))) %>%
        filter(!is.na(Accepted.Description)) %>%
        mutate(across(Compound, ~ paste0("MB_", Compound))) 


# Negative ion mode Compound IDs
negIDs <- neg %>%
        dplyr::select(!matches(c("^X07032024.*\\d\\.1$", "^X07032024.*\\d"))) %>%
        filter(!is.na(Accepted.Description)) %>%
        mutate(across(Compound, ~ paste0("MB_", Compound))) %>% 
        relocate(colnames(posIDs))

# Feature ID table from Pos & Neg ion mode
featIDs <- rbind(posIDs, negIDs)


# Combine Positive and Negative datasets:

# Determine features with lower CV in either Positive or Negative ion mode
common_feats <- intersect(
        unique(posIDs$Accepted.Description %>% na.omit()),
        unique(negIDs$Accepted.Description %>% na.omit())
)



# Remove features with lower CV of QCs

# Read in QC feature lists
posQC <- read.csv("~/Documents/SHIME Project/SHIME/Metabolomics/SHIME_POS_QC_092424.csv", skip = 2, na.string = c(""))
negQC <- read.csv("~/Documents/SHIME Project/SHIME/Metabolomics/SHIME_NEG_QC_092424.csv", skip = 2, na.string = c(""))


# Tidy up

#Positive ion mode
posQCtidy <- posQC %>%
        dplyr::select(!matches("^X07032024.*\\d\\.1$")) %>% 
        filter(!is.na(Accepted.Description) & is.na(treatment)) %>% 
        dplyr::select(Compound, Accepted.Description, starts_with("X07032024")) %>%
        mutate(Compound = paste0("MB_", Compound)) %>%
        dplyr::filter(Accepted.Description %in% common_feats) 

# CV of positive ion mode QCs        
posQC_CVs <- posQCtidy %>% 
        pivot_longer(cols = starts_with("X07032024"), names_to = 'QC', values_to = 'intensity') %>% 
        group_by(Compound, Accepted.Description) %>% 
        summarize(posCV = ( sd(intensity) / mean(intensity) ) * 100 )

#Negative ion mode
negQCtidy <- negQC %>%
        dplyr::select(!matches("^X07032024.*\\d\\.1$")) %>% 
        filter(!is.na(Accepted.Description) & is.na(treatment)) %>% 
        dplyr::select(Compound, Accepted.Description, starts_with("X07032024")) %>%
        mutate(Compound = paste0("MB_", Compound)) %>%
        dplyr::filter(Accepted.Description %in% common_feats)

# CV of negation ion mode QCs
negQC_CVs <- negQCtidy %>% 
        pivot_longer(cols = starts_with("X07032024"), names_to = 'QC', values_to = 'intensity') %>% 
        group_by(Compound, Accepted.Description) %>% 
        summarize(negCV = ( sd(intensity) / mean(intensity) ) * 100 )


QC_CVs <- left_join(posQC_CVs %>% ungroup() %>% dplyr::select(-Compound), 
                    negQC_CVs %>% ungroup() %>% dplyr::select(-Compound)) %>% 
        mutate(to_rm = ifelse(posCV < negCV, 'NEG', 'POS'))



# Wrangle and tidy

# Positive ion mode - includes annotated and unannotated (9023 features, 6968 feature (q<0.05), 106 annotations)
pos_tidy <- pos %>%
        dplyr::select(!matches("^X07032024.*\\d\\.1$")) %>%
        filter(!is.na(anova.q.0.05) & !is.na(Accepted.Description) & is.na(treatment)) %>% 
        dplyr::select(Compound, starts_with("X07032024")) %>% 
        pivot_longer(-Compound) %>% 
        pivot_wider(name, names_from = 'Compound', values_from = 'value') %>% 
        rename(run.id = name) %>% 
        left_join(posids) %>% 
        relocate(sample.id) %>% 
        dplyr::select(-run.id) %>% 
        rename_with(~paste0('MB_', .x), .cols = 2:length(colnames(.))) %>% 
        dplyr::select(-c(
                posQC_CVs %>% filter(Accepted.Description %in%
                                             c(QC_CVs %>% 
                                                       filter(to_rm == 'POS') %>% 
                                                       pull(Accepted.Description))) %>% 
                        pull(Compound)
        ))
# Scale data
pos_scaled <- pos_tidy %>% 
        column_to_rownames('sample.id') %>% 
        mutate(across(everything(), ~ .x + 1)) %>% 
        log_transform() %>% 
        pareto_scale() %>% 
        as.data.frame() %>% 
        rownames_to_column('sample.id')
# Join meta data
pos_full <- right_join(meta, pos_scaled)



# Negative ion mode - includes annotated and unannotated (1755 feature (q<0.05), 206 annotations)
neg_tidy <- neg %>%
        dplyr::select(!matches("^X07032024.*\\d\\.1$")) %>% 
        filter(!is.na(anova.q.0.05) & !is.na(Accepted.Description) & is.na(treatment)) %>%
        dplyr::select(Compound, starts_with("X07032024")) %>% 
        pivot_longer(-Compound,) %>% 
        pivot_wider(name, names_from = 'Compound', values_from = 'value') %>% 
        rename(run.id = name) %>% 
        left_join(negids) %>% 
        relocate(sample.id) %>% 
        dplyr::select(-run.id) %>% 
        rename_with(~paste0('MB_', .x), .cols = 2:length(colnames(.))) %>% 
        dplyr::select(-c(
                negQC_CVs %>% filter(Accepted.Description %in%
                                             c(QC_CVs %>% 
                                                       filter(to_rm == 'NEG') %>% 
                                                       pull(Accepted.Description))) %>% 
                        pull(Compound)
        ))

# Scale data
neg_scaled <- neg_tidy %>% 
        column_to_rownames('sample.id') %>% 
        mutate(across(everything(), ~ .x + 1)) %>% 
        log_transform() %>% 
        pareto_scale() %>% 
        as.data.frame() %>% 
        rownames_to_column('sample.id')

# Join into one big dataset
total_scaled <- left_join(pos_full, neg_scaled)


# Helper Function (reorder data for anova)
anova_helper <- function(dat, dep, indep, condA = 'control', condB = 'treatment'){
        value <- dat[[dep]][sapply(dat[[indep]], function(x) x %in% c(condA, condB))]
        cond <- dat[[indep]][sapply(dat[[indep]], function(x) x %in% c(condA, condB))]
        M <- data.frame(value, cond)
        aov(value ~ cond, data = M)
}

# ANOVA before and after XN treatment
pairwise_by_reactor <- total_scaled %>% 
        column_to_rownames('sample.id') %>% 
        dplyr::select(-date, -day, -chron_order, -to_consider, -week) %>% 
        pivot_longer(starts_with('MB'), names_to = 'Compound', values_to = 'intensity') %>% 
        group_by(SHIME, reactor, Compound) %>% 
        nest() %>%
        mutate(anova = purrr::map(data, function(x) aov(intensity ~ period, data = x))) %>%
        mutate(pval = purrr::map_dbl(anova, function(x) summary(x)[[1]][[5]][[1]])) %>% 
        ungroup()


# Feature Matrix with corrected p-values
sig_feat_mat <- pairwise_by_reactor %>% 
        dplyr::select(SHIME, reactor, Compound, pval) %>% 
        pivot_wider(names_from = c('SHIME', 'reactor'), values_from = pval) %>% 
        left_join(featIDs %>% dplyr::select(Compound, Accepted.Description, compound.class, 
                                            m.z, amino.acid.metabolism, N.acetylation, Adducts, 
                                            Formula, Score, Fragmentation.Score, Mass.Error..ppm., 
                                            Isotope.Similarity), by = 'Compound') %>% 
        filter(!compound.class %in% c('di/tripeptide', 'steroid', 'eicosanoid','other', 'lipid')) %>% 
        mutate(across(A_AC:B_DC, ~replace_na(., 0 ))) 



# Log2 Fold Change for heatmap ------------------------------------

# Join data prior to normalization/transformation
total_prenorm <- left_join(pos_tidy, neg_tidy) %>% left_join(meta) %>% relocate(sample.id, SHIME, reactor, period)


# Log2 FC table of all features
featl2fc <- total_prenorm %>% 
        pivot_longer(starts_with('MB'), names_to = 'Compound', values_to = 'intensity') %>% 
        mutate(intensity = intensity + 1) %>% 
        group_by(SHIME, reactor, Compound) %>% 
        nest() %>% 
        mutate(l2fc = purrr::map(data, 
                                 function(x) 
                                         l2fcmin(x, 'intensity', 'period', 'control', 'treatment'))) %>% 
        dplyr::select(-data) %>% 
        pivot_wider(names_from = c('SHIME', 'reactor'), values_from = 'l2fc') %>% 
        arrange(Compound) %>% 
        mutate(across(everything(), ~unlist(.))) %>% 
        mutate(across(everything(), ~round(., 4))) %>% 
        relocate(A_AC, A_TC, A_DC, B_AC, B_TC, B_DC) %>% 
        left_join(featIDs %>% dplyr::select(Compound, compound.class))


# Figure 5a: Log2 Fold Change of Bile Acid Metabolites -------------------------

# Subset bile acids
bacid_sig_mat <- sig_feat_mat %>% 
        filter(compound.class %in% c('Conjugated primary bile acid', 'Conjugated secondary bile acid', 
                                     'Primary bile acid', 'Secondary bile acid')) %>% 
        mutate(across(A_AC:B_DC, ~replace_na(., 0 ))) %>%
        mutate(across(A_AC:B_DC, ~ p.adjust(.x, method = 'bonferroni'))) %>% 
        arrange(Compound) %>% 
        column_to_rownames('Compound') %>% 
        arrange(compound.class)

# significance symbol map overlay 
bacid_overlay <- bacid_sig_mat %>% 
        dplyr::select(A_AC:B_DC) %>% 
        mutate(across(everything(), ~ifelse(.x <= 0.1 & .x > 0.05, "*",
                                            ifelse(.x <= 0.05 & .x > 0.01, "**", 
                                                   ifelse(.x <= 0.01 & .x > 0, "***", "")))))

# Extract bile acids from log2 fold change
bacid_fc <- featl2fc %>% 
        filter(Compound %in% rownames(bacid_sig_mat)) %>% 
        arrange(compound.class) %>%
        dplyr::select(-compound.class) %>% 
        column_to_rownames('Compound')


# Annotations for columns
col_annotation <- data.frame(SHIME = c(rep('Healthy', time = 3), rep('Crohns', times = 3)),
                             Reactor = rep(c('AC', 'TC', 'DC'), times = 2))
rownames(col_annotation) <- colnames(bacid_fc) 

# Annotations for rows
row_annotation <- bacid_sig_mat %>% 
        dplyr::select(compound.class) %>% 
        rename('Bile acid class' = compound.class)


# Heatmap of Log2 FC 
# [SAVE 1400w x 1400h]
pheatmap::pheatmap(bacid_fc, 
                   color = colorRampPalette(rev(brewer.pal(9, 'RdYlBu')))(24),
                   cluster_rows = F,
                   cluster_cols = F, 
                   display_numbers = bacid_overlay,
                   number_color = 'gray10',
                   fontsize_number = 20,
                   labels_row = bacid_sig_mat$Accepted.Description,
                   annotation_row = row_annotation,
                   border_color = 'grey30',
                   breaks = seq(-8,8, length.out = 24),
                   cellwidth = 35,
                   cellheight = 25,
                   labels_col = c("H-SHIME (AC)", "H-SHIME (TC)", "H-SHIME (DC)",
                                  "D-SHIME (AC)", "D-SHIME (TC)", "D-SHIME (DC)"),
                   fontsize = 16,
                   gaps_row = c(2, 9, 15),
                   angle_col = '315',
                   gaps_col = c(3),
                   legend_breaks = c(-8, -4, 0, 4, 8),
                   legend_labels = c('-8', '-4', '0', '4', '8+'))



# TAUROCHOLIC ACID -------------------------------------------------------------

# metabolomics tidy data
total_tidy <- right_join(meta, pos_tidy) %>% 
        left_join(neg_tidy)


labels <- c('A_AC_control', 'A_TC_control', 'A_DC_control',
            'B_AC_control', 'B_TC_control', 'B_DC_control',
            'A_AC_treatment', 'A_TC_treatment', 'A_DC_treatment',
            'B_AC_treatment', 'B_TC_treatment', 'B_DC_treatment')

# Reorder
tca <- total_tidy %>% 
        dplyr::select(sample.id, SHIME, reactor, date, day, period, week, 'MB_21.66_516.2983m/z') %>% 
        rename(tca = 'MB_21.66_516.2983m/z') %>% 
        mutate(sampreact = paste0(SHIME, "_", reactor)) %>% 
        mutate(cond = paste0(SHIME, '_', reactor, '_', period)) %>% 
        mutate(across(cond, ~ factor(.x, levels = labels))) 



# color palette for SHIME reactors
rct_pal <- c('#CEBBDD' , '#754B95', '#484D6D', '#F3B146',  '#E67B35',  '#EE6352')
names(rct_pal) <- c('A_AC', 'A_TC', 'A_DC', 'B_AC', 'B_TC', 'B_DC')


# Plot TCA intensity pre- and post-XN 
# [SAVE 1000w x 800h]
ggplot(tca, aes(period, tca, color = sampreact))+
        geom_boxplot(alpha = 0.2, size = 1.5, width = 2) +
        #geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.1)) + 
        scale_x_discrete(labels = c('Pre-XN', 'Post-XN')) +
        scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
        scale_fill_manual(values = rct_pal,
                          labels = c('H-SHIME, AC', 'H-SHIME, TC', 'H-SHIME, DC', 
                                     'D-SHIME, AC', 'D-SHIME, TC', 'D-SHIME, DC'),
                          name = 'SHIME, Reactor') + 
        scale_color_manual(values = rct_pal,
                           labels = c('H-SHIME, AC', 'H-SHIME, TC', 'H-SHIME, DC', 
                                      'D-SHIME, AC', 'D-SHIME, TC', 'D-SHIME, DC'),
                           name = 'SHIME, Reactor') + 
        labs(y = 'Intensity', x = '', title = 'Taurocholic acid') +
        theme_bw()+
        theme(axis.title=element_text(size=20),
              axis.text=element_text(size=18),
              strip.text = element_blank(),
              strip.background = element_blank(),
              plot.title = element_text(size=24, hjust = 0.5)) +
        theme(legend.title = element_text(size = 18),
              legend.text = element_text(size = 18)) 



## Figure 5c: BSH Inhibiton ----------------

# BSH inhibition data
bsh <- read.csv("~/Documents/SHIME Project/SHIME/Data/bsh_results_copy.csv") %>% 
        mutate_at('sample', ~ factor(.x, levels = c('ctrl', 'dmso_ctrl', 'XN_1uM', 'XN_10uM', 'XN_50uM'))) 

# Pivot wide
bsh.pivot <- bsh %>% 
        dplyr::select(sample, replicate, CA_d4) %>% 
        pivot_wider(names_from = sample, values_from = CA_d4)

# t-tests compared to DMSO control
t.test(bsh.pivot$dmso_ctrl, bsh.pivot$XN_1uM, data = bsh.pivot) # p=1
t.test(bsh.pivot$dmso_ctrl, bsh.pivot$XN_10uM, data = bsh.pivot) # p=0.05
t.test(bsh.pivot$dmso_ctrl, bsh.pivot$XN_50uM, data = bsh.pivot) # p=0.02

# Summary stats
bsh_summ <- bsh %>% 
        mutate_at('CA_d4', ~ .x/2) %>% 
        group_by(sample) %>% 
        summarize(mCAD4 = mean(CA_d4),
                  SD = sd(CA_d4)) %>% 
        filter(sample != 'ctrl')

# CA-d4 generation by stepwise XN exposure 
# [SAVE 1000w x 800h]
ggplot(bsh_summ, aes(sample, mCAD4)) + 
        geom_bar(stat = 'identity', color = 'gray20', fill = 'white', size = 2, show.legend = F, width = 0.4) +
        geom_errorbar(aes(ymin = mCAD4 - SD, ymax = mCAD4 + SD), size = 2, width = 0.15, color = 'gray20') + 
        labs(y = 'CA-d4 (nM/mg protein per min)', x = "") +
        theme_bw() +
        scale_x_discrete(labels = c('Vehicle', '1 μM', '10 μM', '50 μM')) +
        geom_text(aes(x = 'XN_10uM', y = 1.5, label = '*'), size = 14) +
        geom_text(aes(x = 'XN_50uM', y = 1.5, label = '*'), size = 14) + 
        theme(axis.title=element_text(size=22),
              axis.text=element_text(size=20),
              strip.text = element_blank(),
              strip.background = element_blank())




