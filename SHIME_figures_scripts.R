

# Manuscript for Molecular Nutrition and Food Research Special Issue Polyphenols Application
# Project: SHIME
# November 2024

# Library ----------------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggrepel)
library(ggvenn)
library(phyloseq)
library(vegan)
library(RColorBrewer)


# Data + Environment -----------------------------------------------------------


# Load phyloseq object
shime_ps <- readRDS("~/Documents/16S/SHIME/shime_ps.rds")


# Create placeholder names for un-annotated genera using family
renames <- rownames(tax_table(shime_ps)[is.na(tax_table(shime_ps)[, 'Genus'])])
taxdf <- tax_table(shime_ps)[renames,]
renamed_genus <- unname(sapply(taxa_names(taxdf), function(x) paste0('f_', taxdf[x, 'Family'], '_', x)))
tax_table(shime_ps)[renames, 'Genus'] <- renamed_genus
# Rename NA's from family if present
rename_order <- tax_table(shime_ps) %>% data.frame() %>% filter(., grepl("f_NA", Genus)) %>% rownames()
taxadf <- tax_table(shime_ps)[rename_order,]
renamed_from_order <- unname(sapply(taxa_names(taxadf), function(x) paste0('o_', taxadf[x, 'Order'], '_', x)))
# Rename NA's from Order if present
tax_table(shime_ps)[rename_order, 'Genus'] <- renamed_from_order



# Filter out taxa seen fewer than 10 times in less than 10% of samples
asv_counts <- shime_ps %>% filter_taxa(function(x) sum(x > 10) > (0.1*length(x)), TRUE) # ASVs: 718 to 191
#Filter out low abundance (>1e-5) taxa >>> does not change number of ASVs 
counts_filtered <- asv_counts %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE) # ASVs = 191
# Transform to relative abundance
asv_relab <- asv_counts %>% transform_sample_counts(function(x) x / sum(x))


# Meta Data
meta<- read.csv("~/Documents/SHIME/Data/SHIME_metadata.csv") %>% 
        mutate(across(reactor, ~ factor(.x, levels = c('AC', 'TC', 'DC'))))

# XN Metabolite Data
xn <- read.csv("~/Documents/SHIME/Data/SHIME_xn_metabs.csv") %>% 
        column_to_rownames("sample.id") %>% 
        rename_with(~paste0("XN_", .x)) %>% 
        rownames_to_column("sample.id")

# SCFA Data
scfa <- read.csv("~/Documents/SHIME/Data/SHIME_SCFA_Results_mM.csv") %>% 
        column_to_rownames("sample.id") %>% 
        rename_with(~paste0("FA_", .x)) %>% 
        rownames_to_column("sample.id")



# Functions

# Funciton to calculate standard error
ser <- function(x){
        sd(x)/sqrt(length(x))
}



# Figure 2: XN Metabolites -----------------------------------------------------

# Tidy
xn_combined <- right_join(meta, xn) %>%
        mutate(shime_react = paste0(SHIME, '_', reactor)) %>% 
        pivot_longer(cols = c(starts_with('XN'), -XN_XN), names_to = 'metabolite', values_to = 'conc') %>% 
        mutate(across(metabolite, ~ factor(.x, 
                                           levels = c('XN_IXN', 'XN_DXN', 'XN_x8PN', 'XN_DDXN', 'XN_x6PN'))))


# Calculate mean and st err by week for each metabolite:
xn_reorder <- xn_combined %>% 
        dplyr::filter(week %in% c(5, 6, 7, 8)) %>% 
        group_by(SHIME, reactor, week, metabolite) %>% 
        summarize(mConc = mean(conc),
                  mErr = ser(conc))


xn_pal <- c( '#B6D886', '#FD8B8B', '#F56147', '#F7C059',  '#529163')
names(xn_pal) <- c('XN_IXN', 'XN_DXN', 'XN_x8PN', 'XN_DDXN', 'XN_x6PN')

xn_shimeA <- xn_reorder %>% 
        dplyr::filter(SHIME == 'A')

xn_shimeB <- xn_reorder %>% 
        dplyr::filter(SHIME == 'B')

# Stacked bar chart of XN Metabolties (SHIME A) [300w x 1000h]
ggplot(xn_shimeA, aes(week, mConc, fill = metabolite, group = SHIME)) +
        geom_bar(stat = 'identity') + 
        facet_wrap(~ reactor, ncol = 1) +
        scale_fill_manual(values = xn_pal, 
                          labels = c('IXN', 'DXN', '8PN', 'DDXN', '6PN'), 
                          name = 'Metabolite',
                          guide = 'none') +
        theme_bw() +
        scale_y_continuous(limits = c(0, 150)) + 
        labs(x = 'week', y = '') +
        theme(axis.title=element_text(size=24),
              axis.text=element_text(size=22),
              strip.text = element_text(size=24))

# Stacked bar chart of XN Metabolties (SHIME B) [300w x 1000h]
ggplot(xn_shimeB, aes(week, mConc, fill = metabolite, group = SHIME)) +
        geom_bar(stat = 'identity') + 
        facet_wrap(~ reactor, ncol = 1) +
        scale_fill_manual(values = xn_pal, 
                          labels = c('IXN', 'DXN', '8PN', 'DDXN', '6PN'), 
                          name = 'Metabolite',
                          guide = 'none') +
        theme_bw() +
        scale_y_continuous(limits = c(0, 150)) + 
        labs(x = 'week', y = '') +
        theme(axis.title=element_text(size=24),
              axis.text=element_text(size=22),
              strip.text = element_text(size=24))


# Longitudinal plots

shimeA <- xn_combined %>% 
        filter(SHIME == 'A')

shimeB <- xn_combined %>% 
        filter(SHIME == 'B')

# shime A [600 x 1000]
ggplot(shimeA, aes(day, conc, fill = metabolite, group = metabolite, color = metabolite)) + 
        geom_path(size = 1.5) + 
        facet_wrap(~reactor, ncol = 1) + 
        scale_color_manual(values = xn_pal, 
                           labels = c('IXN', 'DXN', '8PN', 'DDXN', '6PN'), 
                           name = 'Metabolite',
                           guide = 'none') +
        scale_fill_manual(values = xn_pal) +
        theme_bw() +
        scale_y_continuous(limits = c(0, 250)) + 
        labs(x = 'day', y = 'nM') +
        theme(axis.title=element_text(size=24),
              axis.text=element_text(size=22),
              strip.text = element_text(size=24))

# shime B [600 x 1000]
ggplot(shimeB, aes(day, conc, fill = metabolite, group = metabolite, color = metabolite)) + 
        geom_path(size = 1.5) + 
        facet_wrap(~reactor, ncol = 1) + 
        scale_color_manual(values = xn_pal, 
                           labels = c('IXN', 'DXN', '8PN', 'DDXN', '6PN'), 
                           name = 'Metabolite',
                           guide = 'none') +
        scale_fill_manual(values = xn_pal) +
        theme_bw() +
        scale_y_continuous(limits = c(0, 250)) + 
        labs(x = 'day', y = 'nM') +
        theme(axis.title=element_text(size=24),
              axis.text=element_text(size=22),
              strip.text = element_text(size=24))


L <- ggplot(xn_shimeA, aes(week, mConc, fill = metabolite, group = SHIME)) +
        geom_bar(stat = 'identity') + 
        facet_wrap(~ reactor, ncol = 1) +
        scale_fill_manual(values = xn_pal, 
                          labels = c('IXN', 'DXN', '8PN', 'DDXN', '6PN'), 
                          name = 'Metabolite') +
        theme_bw() +
        scale_y_continuous(limits = c(0, 150)) + 
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

# [700w x 600h]
as_ggplot(leg)



# Figure 3: SCFA concentrations ------------------------------------------------

fa_total <- left_join(meta, scfa) %>% 
        column_to_rownames("sample.id") %>% 
        mutate(shime_react = paste0(SHIME, "_", reactor)) %>% 
        dplyr::filter(period != 'stabilization')


rct_pal <- c('#CEBBDD' , '#754B95', '#484D6D', '#F3B146',  '#E67B35',  '#EE6352')
names(rct_pal) <- c('A_AC', 'A_TC', 'A_DC', 'B_AC', 'B_TC', 'B_DC')


fa1 <- ggplot(fa_total, aes(day, FA_Total, group = shime_react, color = shime_react)) +
        geom_path(size = 1.8) + 
        geom_vline(xintercept = 30, linetype = 2, size = 1) +
        annotate('text', x = 33, y = 0, label = 'treatment', size = 8) +
        theme_bw() +
        scale_y_continuous(limits = c(0, 60)) + 
        scale_x_continuous(breaks = c(15, 23, 29, 36, 43, 50), labels = c(0, 7, 14, 21, 28, 35)) + 
        labs(x = 'day', y = 'Total SCFA (mM)') +
        theme(axis.title=element_text(size=28),
              axis.text=element_text(size=24),
              strip.text = element_text(size=26)) +
        scale_color_manual(values = rct_pal,
                           labels = c('H-SHIME, AC', 'H-SHIME, TC', 'H-SHIME, DC', 
                                      'C-SHIME, AC', 'C-SHIME, TC', 'C-SHIME, DC'),
                           name = 'SHIME, Reactor',
                           guide = 'none') 


fa2 <- ggplot(fa_total, aes(day, FA_Acetic, group = shime_react, color = shime_react)) +
        geom_path(size = 1.8) + 
        geom_vline(xintercept = 30, linetype = 2, size = 1) +
        annotate('text', x = 33, y = 0, label = 'treatment', size = 8) +
        theme_bw() +
        scale_y_continuous(limits = c(0, 40)) + 
        scale_x_continuous(breaks = c(15, 23, 29, 36, 43, 50), labels = c(0, 7, 14, 21, 28, 35)) + 
        labs(x = 'day', y = 'Acetic acid (mM)') +
        theme(axis.title=element_text(size=28),
              axis.text=element_text(size=24),
              strip.text = element_text(size=26)) +
        scale_color_manual(values = rct_pal,
                           labels = c('H-SHIME, AC', 'H-SHIME, TC', 'H-SHIME, DC', 
                                      'C-SHIME, AC', 'C-SHIME, TC', 'C-SHIME, DC'),
                           name = 'SHIME, Reactor',
                           guide = 'none') 

fa3 <- ggplot(fa_total, aes(day, FA_Propionic, group = shime_react, color = shime_react)) +
        geom_path(size = 1.8) + 
        geom_vline(xintercept = 30, linetype = 2, size = 1) +
        annotate('text', x = 33, y = 0, label = 'treatment', size = 8) +
        theme_bw() +
        scale_y_continuous(limits = c(0, 35)) + 
        scale_x_continuous(breaks = c(15, 23, 29, 36, 43, 50), labels = c(0, 7, 14, 21, 28, 35)) + 
        labs(x = 'day', y = 'Propionic acid (mM)') +
        theme(axis.title=element_text(size=28),
              axis.text=element_text(size=24),
              strip.text = element_text(size=26)) +
        scale_color_manual(values = rct_pal,
                           labels = c('H-SHIME, AC', 'H-SHIME, TC', 'H-SHIME, DC', 
                                      'C-SHIME, AC', 'C-SHIME, TC', 'C-SHIME, DC'),
                           name = 'SHIME, Reactor',
                           guide = 'none') 

fa4 <- ggplot(fa_total, aes(day, FA_Butyric, group = shime_react, color = shime_react)) +
        geom_path(size = 1.8) + 
        geom_vline(xintercept = 30, linetype = 2, size = 1) +
        annotate('text', x = 33, y = 0, label = 'treatment', size = 8) +
        theme_bw() +
        scale_y_continuous(limits = c(0, 35)) + 
        scale_x_continuous(breaks = c(15, 23, 29, 36, 43, 50), labels = c(0, 7, 14, 21, 28, 35)) + 
        labs(x = 'day', y = 'Butyric acid (mM)') +
        theme(axis.title=element_text(size=28),
              axis.text=element_text(size=24),
              strip.text = element_text(size=26)) +
        scale_color_manual(values = rct_pal,
                           labels = c('H-SHIME, AC', 'H-SHIME, TC', 'H-SHIME, DC', 
                                      'C-SHIME, AC', 'C-SHIME, TC', 'C-SHIME, DC'),
                           name = 'SHIME, Reactor',
                           guide = 'none') 

figure <- ggarrange(fa1, fa2, fa3, fa4,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2,
                    font.label = list(size = 28, color = "black"),
                    common.legend = T)
# [2200 x 1200]
figure


# Extract Legend
L <- ggplot(fa_total, aes(day, FA_Butyric, group = shime_react, color = shime_react)) +
        geom_path(size = 5) + 
        geom_vline(xintercept = 30, linetype = 2, size = 1) +
        annotate('text', x = 33, y = 0, label = 'treatment', size = 8) +
        theme_bw() +
        scale_y_continuous(limits = c(0, 35)) + 
        labs(x = 'day', y = 'Butyric acid (mM)') +
        theme(axis.title=element_text(size=26),
              axis.text=element_text(size=22),
              strip.text = element_text(size=24)) +
        scale_color_manual(values = rct_pal,
                           labels = c('AC', 'TC', 'DC', 
                                      'AC', 'TC', 'DC'),
                           name = '') +
        theme(legend.key.size = unit(2, 'cm'), #change legend key size
              legend.key.height = unit(2, 'cm'), #change legend key height
              legend.key.width = unit(2, 'cm'), #change legend key width
              legend.title = element_text(size=30), #change legend title font size
              legend.text = element_text(size=26), #change legend text font size
              legend.position = 'bottom') + 
        guides(colour = guide_legend(nrow = 1))


leg <- get_legend(L)

# [700w x 600h]
as_ggplot(leg)





# Figure 4: 16S MICROBIOME -----------------------------------------------------



# Figure 4a: Alpha Diversity ---------------------------------------------------

# Subset treatment and control periods
shime_treat <- shime_ps %>% subset_samples(week != '5' & period %in% c('control', 'treatment'))

# Alpha Diversity
meta <- shime_treat %>% sample_data() %>% data.frame() %>% rownames_to_column('sample.id')

adiv <- shime_treat %>%
        estimate_richness(measures = c('Observed', 'Shannon', 'Simpson')) %>% 
        rownames_to_column('sample.id') %>% 
        left_join(., meta)


rct_pal <- c('#CEBBDD' , '#754B95', '#484D6D', '#F3B146',  '#E67B35',  '#EE6352')
names(rct_pal) <- c('A_AC', 'A_TC', 'A_DC', 'B_AC', 'B_TC', 'B_DC')


ad1 <- ggplot(adiv, aes(day, Observed, group = sampreact, color = sampreact)) +
        geom_path(size = 1.8) + 
        geom_vline(xintercept = 30, linetype = 2, size = 1) +
        annotate('text', x = 34, y = 4, label = 'treatment', size = 8) +
        #theme_bw() +
        scale_y_continuous(limits = c(0, 175)) + 
        scale_x_continuous(breaks = c(15, 23, 29, 36, 43, 50), labels = c(0, 7, 14, 21, 28, 35)) + 
        labs(x = '', y = 'Observed') +
        theme(axis.title=element_text(size=26),
              axis.text=element_text(size=24),
              strip.text = element_blank(),
              strip.background = element_blank()) +
        scale_color_manual(values = rct_pal,
                           labels = c('H-SHIME, AC', 'H-SHIME, TC', 'H-SHIME, DC', 
                                      'C-SHIME, AC', 'C-SHIME, TC', 'C-SHIME, DC'),
                           name = 'SHIME, Reactor',
                           guide = 'none') +
        facet_wrap(~SHIME, ncol = 2)

ad2 <- ggplot(adiv, aes(day, Shannon, group = sampreact, color = sampreact)) +
        geom_path(size = 1.8) + 
        geom_vline(xintercept = 30, linetype = 2, size = 1) +
        annotate('text', x = 34, y = 0.4, label = 'treatment', size = 8) +
        #theme_bw() +
        scale_y_continuous(limits = c(0, 4.2)) + 
        scale_x_continuous(breaks = c(15, 23, 29, 36, 43, 50), labels = c(0, 7, 14, 21, 28, 35)) + 
        labs(x = '', y = 'Shannon') +
        theme(axis.title=element_text(size=26),
              axis.text=element_text(size=24),
              strip.text = element_blank(),
              strip.background = element_blank()) +
        scale_color_manual(values = rct_pal,
                           labels = c('H-SHIME, AC', 'H-SHIME, TC', 'H-SHIME, DC', 
                                      'C-SHIME, AC', 'C-SHIME, TC', 'C-SHIME, DC'),
                           name = 'SHIME, Reactor',
                           guide = 'none') +
        facet_wrap(~SHIME, ncol = 2)

ad3 <- ggplot(adiv, aes(day, Simpson, group = sampreact, color = sampreact)) +
        geom_path(size = 1.8) + 
        geom_vline(xintercept = 30, linetype = 2, size = 1) +
        annotate('text', x = 34, y = 0.1, label = 'treatment', size = 8) +
        #theme_bw() +
        scale_y_continuous(limits = c(0, 1.2)) + 
        scale_x_continuous(breaks = c(15, 23, 29, 36, 43, 50), labels = c(0, 7, 14, 21, 28, 35)) + 
        labs(x = 'day', y = 'Simpson') +
        theme(axis.title=element_text(size=26),
              axis.text=element_text(size=24),
              strip.text = element_blank(),
              strip.background = element_blank()) +
        scale_color_manual(values = rct_pal,
                           labels = c('H-SHIME, AC', 'H-SHIME, TC', 'H-SHIME, DC', 
                                      'C-SHIME, AC', 'C-SHIME, TC', 'C-SHIME, DC'),
                           name = 'SHIME, Reactor',
                           guide = 'none') +
        facet_wrap(~SHIME, ncol = 2)

figure <- ggarrange(ad1, ad2, ad3,
                    ncol = 1, nrow = 3,
                    align = 'v',
                    font.label = list(size = 24, color = "black"))
# [1600 x 1000]
figure





# Figure 4b: Aitchison PCA + BiPLOT --------------------------------------------

# Subset treatment and control periods
treat_counts <- asv_counts %>% subset_samples(week != '5' & period %in% c('control', 'treatment'))

clr <- treat_counts %>% 
        rarefy_even_depth(rngseed = 44) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1)

meta <- treat_counts %>% 
        sample_data() %>% 
        data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        mutate(treatment = ifelse(period == 'control', 'pre XN', 'post XN'))

euclid_dist <- clr %>% 
        stats::prcomp() 

comps <- euclid_dist$x %>% 
        as.data.frame %>% 
        rownames_to_column('sample.id') %>% 
        right_join(meta) %>% 
        mutate(across(sampreact, ~ factor(.x, levels = c('A_AC', 'A_TC', 'A_DC', 'B_AC', 'B_TC', 'B_DC'))))

aitch_pca <- ggplot(comps, aes(x = PC1, y = PC2, shape = period, color = sampreact))+
        geom_point(size=8) + 
        labs(x = 'PC1 [46.3%]', y = 'PC2 [11.6%]', shape = '') +
        scale_color_manual(values = rct_pal,
                           labels = c('H-SHIME, AC', 'H-SHIME, TC', 'H-SHIME, DC', 
                                      'C-SHIME, AC', 'C-SHIME, TC', 'C-SHIME, DC'),
                           name = '') +
        theme(axis.title=element_text(size=26),
              axis.text=element_text(size=22),
              strip.text = element_text(size=24)) +
        theme(legend.key.size = unit(1, 'cm'), #change legend key size
              legend.key.height = unit(1, 'cm'), #change legend key height
              legend.key.width = unit(1, 'cm'), #change legend key width
              legend.title = element_text(size=28), #change legend title font size
              legend.text = element_text(size=26),
              legend.position = 'bottom') + 
        guides(colour = guide_legend(nrow = 2))

# BiPLot (add vectors with strongest correlation to PC1 and PC2)

tax <- tax_table(treat_counts) %>% 
        data.frame %>% 
        rownames_to_column('OTU')

melt_combine <- treat_counts %>% 
        psmelt() %>% 
        rename(sample.id = Sample) %>% 
        dplyr::select(OTU, sample.id, Abundance) %>% 
        left_join(comps %>% dplyr::select(sample.id, PC1, PC2), by = 'sample.id')

melt.cor <- melt_combine %>% 
        group_by(OTU) %>% 
        nest() %>% 
        mutate(corX = purrr::map(data, ~ broom::tidy(cor.test(.x$Abundance, .x$PC1, method = 'spearman', exact = F)))) %>%
        mutate(corY = purrr::map(data, ~ broom::tidy(cor.test(.x$Abundance, .x$PC2, method = 'spearman', exact = F)))) %>% 
        dplyr::select(-data) %>% 
        unnest(corX, corY) %>% 
        dplyr::select(OTU, estimate, p.value, estimate1, p.value1)

high_corr <- melt.cor %>% 
        filter(abs(estimate) > 0.84 | abs(estimate1) > 0.70) %>% 
        mutate(estX = estimate * 10) %>% 
        mutate(estY = estimate1 * 10) %>% 
        left_join(tax %>% dplyr::select(OTU, Genus, Species), by = 'OTU') %>% 
        mutate(names = ifelse(is.na(Species), paste0(Genus, "_", OTU), paste(Genus, "_", Species)))

library(ggrepel)

# [1800 x 1200]
aitch_biplot <- aitch_pca +
        geom_segment(data = high_corr, aes(x = 0, xend = estX, y = 0, yend = estY),
                     inherit.aes = F, arrow=arrow(length=unit(0.4,"cm")), alpha=0.75, color="red") + 
        geom_text_repel(data = high_corr, aes(x = estX, y = estY, label = names),
                        inherit.aes = F, alpha = 0.75, color = "red", size = 8)




# Figure 4c: Venn Diagram ------------------------------------------------------

asv_melt <- treat_counts %>% 
        psmelt() %>% 
        dplyr::select(OTU, Sample, Abundance, SHIME, reactor, period) 

A_pre <- asv_melt %>% 
        filter(SHIME == 'A' & period == 'control') %>% 
        filter(Abundance != 0) %>% 
        pull(OTU) %>% 
        unique()

A_post <- asv_melt %>% 
        filter(SHIME == 'A' & period == 'treatment') %>% 
        filter(Abundance != 0) %>%  
        pull(OTU) %>% 
        unique()

B_pre <- asv_melt %>% 
        filter(SHIME == 'B' & period == 'control') %>% 
        filter(Abundance != 0) %>%  
        pull(OTU) %>% 
        unique()

B_post <- asv_melt %>% 
        filter(SHIME == 'B' & period == 'treatment') %>% 
        filter(Abundance != 0) %>%  
        pull(OTU) %>% 
        unique()


venn <- list(
        "H-SHIME (ctrl)" = A_pre,
        "H-SHIME (XN)" = A_post,
        "C-SHIME (XN)" = B_post,
        "C-SHIME (ctrl)" = B_pre
)

library(ggvenn)
ggvenn(
        venn, 
        fill_color = c("#0073C2FF", '#F7C059', '#B6D886', '#F56147'),
        stroke_size = 0.8,
        stroke_linetype = "dashed",
        set_name_size = 4.5,
        show_percentage = F
) 



# Figure 4e: Relative Abundance stacked bar charts -----------------------------

# Family Level
treat_counts <- asv_counts %>% subset_samples(period %in% c('control', 'treatment'))

# First average counts by week, then convert to relative abundance 
fam_ps <- treat_counts %>% 
        tax_glom('Family') %>% 
        psmelt() %>% 
        group_by(SHIME, reactor, week, Family) %>% 
        summarize(mAbund = mean(Abundance)) %>% 
        ungroup(Family) %>% 
        mutate(relabund = mAbund / sum(mAbund)) %>% 
        dplyr::select(-mAbund)


colorcount <- length(unique(fam_ps$Family))
getPalette <- colorRampPalette(brewer.pal(10, 'Set3'))


colorcount <- length(unique(fam_ps$Family))
getPalette <- colorRampPalette(brewer.pal(12, 'Paired'))


# 1400 x 1000
ggplot(fam_ps, aes(x = week, y = relabund, fill = Family, group = SHIME)) +
        geom_bar(stat = 'identity', position = position_fill()) +
        cowplot::theme_cowplot() +
        scale_fill_manual(values = getPalette(colorcount)) +
        facet_grid(reactor ~ SHIME, labeller = labeller(SHIME = c('A' = 'H-SHIME',
                                                                  'B' = 'C-SHIME'))) + 
        guides(fill=guide_legend(ncol=1)) +
        labs(y = 'Relative Abundance') +
        scale_y_continuous(breaks = c(0.00, 1.00), labels = c('0%', '100%')) +
        theme(axis.title=element_text(size=26),
              axis.text=element_text(size=22),
              strip.text = element_text(size=24)) +
        theme(
                legend.title = element_text(size=22), #change legend title font size
                legend.text = element_text(size=18)) #change legend text font size





# Figure 4d: GLMNB Differential Abundance Analysis -----------------------------

# from analysis script
glm.fil <- glmDA %>% 
        dplyr::select(OTU, SHIME, summ, padj) %>% 
        mutate(trt = purrr::map_dbl(summ, function(x) x$coefficients[[2,1]])) %>% 
        filter(padj <= 0.05 & padj > 0) %>% 
        left_join(tax, by = 'OTU')


# For Heatmap 
counts_filtered <- asv_counts %>% subset_samples(week != '5' & period %in% c('control', 'treatment'))

tax <- tax_table(counts_filtered) %>% data.frame() %>% rownames_to_column('OTU') %>% 
        filter(OTU %in% glm.fil$OTU)

# Use counts
asvcnts <- counts_filtered %>% 
        rarefy_even_depth(rngseed = 23) %>% 
        otu_table() %>% 
        data.frame()%>% 
        rownames_to_column('sample.id')
asvmd <- counts_filtered %>% 
        sample_data() %>% 
        data.frame() %>% 
        rownames_to_column('sample.id')
asvcntmd <- left_join(asvmd, asvcnts)

# Use count data for log2FC
asvl2fc <- asvcntmd %>% 
        pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = 'counts') %>% 
        mutate(counts = counts + 1) %>% 
        group_by(SHIME, reactor, ASV) %>% 
        nest() %>% 
        mutate(l2fc = purrr::map(data, function(x) l2fcmin(x, 'counts', 'phase', 'pre', 'post'))) %>% 
        dplyr::select(-data) %>% 
        pivot_wider(names_from = c('SHIME', 'reactor'), values_from = 'l2fc') %>% 
        # Filter significant ASV from GLMNB model
        dplyr::filter(ASV %in% all_of(glm.fil$OTU)) %>% 
        arrange(ASV) %>% 
        column_to_rownames('ASV') %>% 
        mutate(across(everything(), ~as.numeric(.))) %>% 
        mutate(across(everything(), ~round(., 3))) 

#Matrix of sig values
pvals <- glm.fil %>% 
        dplyr::select(OTU, SHIME, padj) %>% 
        pivot_wider(names_from = c('SHIME'), values_from = 'padj') %>% 
        rename(ASV = OTU, B_AC = B, A_AC = A) %>%
        mutate(B_TC = B_AC, B_DC = B_AC, A_TC = A_AC, A_DC = A_AC) %>% 
        relocate(A_AC, A_TC, A_DC, B_AC, B_TC, B_DC) %>% 
        arrange(ASV) %>% 
        column_to_rownames('ASV') %>% 
        mutate(across(everything(), ~replace_na(., 0))) 


#Sig matrix convert to * symbol 
sigmat <- pvals %>% 
        mutate(across(everything(), ~ifelse(.x <= 0.05 & .x > 0.01, "*", 
                                            ifelse(.x <= 0.01 & .x > 0.001, "**",
                                                   ifelse(.x <= 0.001 & .x > 0, "***", "")))))

labels = tax %>% 
        filter(OTU %in% rownames(asvl2fc)) %>% 
        arrange(OTU) %>% 
        mutate(names = ifelse(is.na(Species), paste0(Genus, "_", OTU), paste0(Genus, "_", Species, "_", OTU))) %>% 
        column_to_rownames('OTU')


#Heatmap of log2FC [1000 x 1600]
pheatmap::pheatmap(asvl2fc,
                   color = colorRampPalette(rev(brewer.pal(10, 'RdBu')))(24),
                   cluster_rows = T,
                   cluster_cols = F, 
                   display_numbers = sigmat,
                   number_color = 'grey100',
                   fontsize_number = 16,
                   labels_row = labels$names,
                   annotation_row = labels[, 5, drop = F],
                   breaks = seq(-4,4, length.out = 24),
                   cellwidth = 45,
                   cellheight = 25,
                   labels_col = c("H-SHIME (AC)", "H-SHIME (TC)", "H-SHIME (DC)",
                                  "C-SHIME (AC)", "C-SHIME (TC)", "C-SHIME (DC)"),
                   angle_col = '315',
                   gaps_col = c(3))





# Figure 5: Metabolomics -------------------------------------------------------


# Figure 5a: Principal Component Analysis --------------------------------------

# variable from analysis script - conduct PCA
MBdist <- total_scaled %>% 
        dplyr::select(sample.id, starts_with('MB')) %>% 
        column_to_rownames('sample.id') %>% 
        stats::prcomp()

# extract principal components and join data
MBcomps <- MBdist$x %>% 
        as.data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        right_join(total_scaled %>% dplyr::select(sample.id, SHIME, reactor, period, day)) %>% 
        mutate(shime_react = paste0(SHIME, '_', reactor)) %>% 
        mutate(across(shime_react, ~ factor(.x, level = c('A_AC', 'A_TC', 'A_DC', 'B_AC', 'B_TC', 'B_DC'))))

# Plot PCA
MB_pca <- ggplot(MBcomps, aes(PC1, PC2, shape = period, color = shime_react)) +
        geom_point(size=8) + 
        labs(x = 'PC1 [32.2%]', y = 'PC2 [22.0%]', shape = '') +
        scale_color_manual(values = rct_pal,
                           labels = c('H-SHIME, AC', 'H-SHIME, TC', 'H-SHIME, DC', 
                                      'C-SHIME, AC', 'C-SHIME, TC', 'C-SHIME, DC'),
                           name = 'SHIME, Reactor') +
        theme(axis.title=element_text(size=26),
              axis.text=element_text(size=22),
              strip.text = element_text(size=24)) +
        theme(legend.key.size = unit(2, 'cm'), #change legend key size
              legend.key.height = unit(2, 'cm'), #change legend key height
              legend.key.width = unit(2, 'cm'), #change legend key width
              legend.title = element_text(size=30), #change legend title font size
              legend.text = element_text(size=26)) #change legend text font size





# Figure 5b: ChemRICH Analysis Plot --------------------------------------------

# Read in ChemRICH results
CHR <- read.csv("~/Documents/SHIME/Metabolomics/CHRplot.csv") %>% 
        dplyr::filter(Cluster.name != 'dipeptides') 


# [1200 x 1000]
ggplot(CHR, aes(x = xlogp, y = -log(p.values))) +
        geom_point(aes(size = Cluster.size, color = Increased.ratio)) +
        scale_color_gradient(low = '#194CA5', high = 'red', limits = c(0,1)) +
        scale_size(range = c(10, 40)) +
        scale_x_continuous(" Median XlogP (lipophilicity) ") +
        theme_bw() +
        geom_label_repel(aes(label = Cluster.name), color = "gray20",family="Arial",data=subset(CHR, Cluster.size>2),
                         force = 12, size = 7, force_pull = 4)+
        #theme(text=element_text(family="Arial Black")) +
        theme(
                plot.title = element_text(size=30, hjust = 0.5),
                axis.title.x = element_text(size=27),
                axis.title.y = element_text(size=27, angle=90),
                panel.grid.major = element_blank(), # switch off major gridlines
                panel.grid.minor = element_blank(), # switch off minor gridlines
                legend.position = "none", # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
                legend.title = element_blank(), # switch off the legend title
                legend.text = element_text(size=12),
                legend.key.size = unit(1.5, "lines"),
                legend.key = element_blank(), # switch off the rectangle around symbols in the legend
                legend.spacing = unit(.05, "cm"),
                axis.text.x = element_text(size=20,angle = 0, hjust = 1),
                axis.text.y = element_text(size=25,angle = 0, hjust = 1)
        )



# Figure 5c: Metabolite log2 FC Heatmaps ---------------------------------------


# pairwise ttest helper function
ttest_helper <- function(dat, condition, measure, condA = 'control', condB = 'treatment'){
        cA <- dat[[measure]][sapply(dat[[condition]], function(x) x == condA)]
        cB <- dat[[measure]][sapply(dat[[condition]], function(x) x == condB)]
        M <- data.frame(ctrl = cA, trt = cB)
        t.test(Pair(trt, ctrl) ~ 1, data = M)
        
}


# Pairwise ttest before and after XN treatment
pairwise_by_reactor <- total_scaled %>% 
        pivot_longer(starts_with('MB'), names_to = 'Compound', values_to = 'intensity') %>% 
        group_by(SHIME, reactor, Compound) %>% 
        nest() %>% 
        mutate(test = purrr::map(data, function(x) ttest_helper(dat=x, 'period', 'intensity'))) %>% 
        #mutate(test = purrr::map(data, function(x) t.test(x$intensity ~ x$period, alternative = 'two.sided', paired = T))) %>% 
        mutate(pval = purrr::map_dbl(test, function(x) x$p.value)) %>% 
        ungroup()


# Feature Matrix with corrected p-values
sig_feat_mat <- pairwise_by_reactor %>% 
        dplyr::select(SHIME, reactor, Compound, pval) %>% 
        pivot_wider(names_from = c('SHIME', 'reactor'), values_from = pval) %>% 
        left_join(featIDs %>% dplyr::select(Compound, Accepted.Description, compound.class, 
                                            m.z, amino.acid.metabolism, N.acetylation, Adducts, 
                                            Formula, Score, Fragmentation.Score, Mass.Error..ppm., 
                                            Isotope.Similarity), by = 'Compound') %>% 
        filter(compound.class != "di/tripeptide") %>% 
        filter(!compound.class %in% c('steroid', 'eicosanoid','other'))



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



# Subset Amino Acids (supplemental figure S1) ----------------------------------

# Subset Amino Acids
amino_sig_mat <- sig_feat_mat %>% 
        filter(!is.na(amino.acid.metabolism)) %>% 
        filter(is.na(N.acetylation)) %>% 
        mutate(across(A_AC:B_DC, ~replace_na(., 0 ))) %>%
        mutate(across(A_AC:B_DC, ~ p.adjust(.x, method = 'BH'))) %>% 
        arrange(Compound) %>% 
        column_to_rownames('Compound') %>% 
        arrange(compound.class)


# Matrix to overlay significant symbols (amino acids)
amino_overlay <- amino_sig_mat %>% 
        dplyr::select(A_AC:B_DC) %>% 
        mutate(across(everything(), ~ifelse(.x <= 0.1 & .x > 0.05, "*",
                                            ifelse(.x <= 0.05 & .x > 0.01, "**", 
                                                   ifelse(.x <= 0.01 & .x > 0, "***", "")))))
amino_fc <- featl2fc %>% 
        filter(Compound %in% rownames(amino_sig_mat)) %>% 
        arrange(compound.class) %>%
        dplyr::select(-compound.class) %>% 
        column_to_rownames('Compound')


# Annotations for columns
col_annotation <- data.frame(SHIME = c(rep('Healthy', time = 3), rep('Crohns', times = 3)),
                             Reactor = rep(c('AC', 'TC', 'DC'), times = 2))
rownames(col_annotation) <- colnames(amino_fc) 

# Annotations for rows
row_annotation <- amino_sig_mat %>% 
        dplyr::select(compound.class) %>% 
        rename('Amino acid precursor' = compound.class)

# Heatmap of Log2 FC [1600 x 2000]
pheatmap::pheatmap(amino_fc, 
                   color = colorRampPalette(rev(brewer.pal(9, 'RdYlBu')))(24),
                   cluster_rows = F,
                   cluster_cols = F, 
                   display_numbers = amino_overlay,
                   number_color = 'gray10',
                   fontsize_number = 20,
                   labels_row = amino_sig_mat$Accepted.Description,
                   annotation_row = row_annotation,
                   border_color = 'grey30',
                   breaks = seq(-5.5,5.5, length.out = 24),
                   cellwidth = 35,
                   cellheight = 20,
                   labels_col = c("H-SHIME (AC)", "H-SHIME (TC)", "H-SHIME (DC)",
                                  "C-SHIME (AC)", "C-SHIME (TC)", "C-SHIME (DC)"),
                   fontsize = 16,
                   gaps_row = c(14, 17, 28, 37, 38, 55),
                   angle_col = '315',
                   gaps_col = c(3))




# Figure 5c: Log2 Fold Change of Bile Acid Metabolites -------------------------


# Subset amino acids
bacid_sig_mat <- sig_feat_mat %>% 
        filter(compound.class %in% c('Conjugated primary bile acid', 'Conjugated secondary bile acid', 
                                     'Primary bile acid', 'Secondary bile acid')) %>% 
        mutate(across(A_AC:B_DC, ~replace_na(., 0 ))) %>%
        mutate(across(A_AC:B_DC, ~ p.adjust(.x, method = 'BH'))) %>% 
        arrange(Compound) %>% 
        column_to_rownames('Compound') %>% 
        arrange(compound.class)


bacid_overlay <- bacid_sig_mat %>% 
        dplyr::select(A_AC:B_DC) %>% 
        mutate(across(everything(), ~ifelse(.x <= 0.1 & .x > 0.05, "*",
                                            ifelse(.x <= 0.05 & .x > 0.01, "**", 
                                                   ifelse(.x <= 0.01 & .x > 0, "***", "")))))


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


# Heatmap of Log2 FC [1400 x 1400]
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
                   breaks = seq(-6,6, length.out = 24),
                   cellwidth = 35,
                   cellheight = 25,
                   labels_col = c("H-SHIME (AC)", "H-SHIME (TC)", "H-SHIME (DC)",
                                  "C-SHIME (AC)", "C-SHIME (TC)", "C-SHIME (DC)"),
                   fontsize = 16,
                   gaps_row = c(2, 9, 15),
                   angle_col = '315',
                   gaps_col = c(3))




# Figure 6: BSH Activity -------------------------------------------------------


# Load results from BSH activity assay
bsh <- read.csv("~/Documents/SHIME/Data/bsh_results.csv") %>% 
        mutate_at('sample', ~ factor(.x, levels = c('ctrl', 'dmso_ctrl', 'XN_1uM', 'XN_10uM', 'XN_50uM')))

bsh_summ <- bsh %>% 
        mutate_at('CA_d4', ~ .x/2) %>% 
        group_by(sample) %>% 
        summarize(mCAD4 = mean(CA_d4),
                  SD = sd(CA_d4))


con.pal <- c('#96CCF5', '#F3856F', '#B1DCBA' , '#B1DCBA' , '#B1DCBA' )
names(con.pal) <- c('ctrl', 'dmso_ctrl', 'XN_1uM', 'XN_10uM', 'XN_50uM')


# 1000 x 800
ggplot(bsh_summ, aes(sample, mCAD4, fill = sample)) + 
        geom_bar(stat = 'identity', color = 'black', size = 1.6, show.legend = F, width = 0.5)+
        geom_errorbar(aes(ymin = mCAD4 - SD, ymax = mCAD4 + SD), size = 1.6, width = 0.2) + 
        labs(y = 'CA-d4 (nM/mg protein per min)', x = "") +
        theme_bw() +
        scale_x_discrete(labels = c('Control', 'Vehicle', 'XN (1 μM)', 'XN (10 μM)', 'XN (50 μM)')) +
        scale_fill_manual(values = con.pal) +
        #geom_segment(aes(x = 'dmso_ctrl', xend = 'XN_10uM', y = 2, yend = 2), size = 1.5) +
        #geom_segment(aes(x = 'ctrl', xend = 'XN_10uM', y = 2.2, yend = 2.2), size = 1.5) +
        #geom_segment(aes(x = 'ctrl', xend = 'XN_50uM', y = 2.4, yend = 2.4), size = 1.5) +
        geom_text(aes(x = 'XN_10uM', y = 1.5, label = '*'), size = 22) +
        geom_text(aes(x = 'XN_50uM', y = 1.5, label = '*'), size = 22) + 
        theme(axis.title=element_text(size=32),
              axis.text=element_text(size=28),
              strip.text = element_blank(),
              strip.background = element_blank(),
              axis.line = element_line(colour = 'black', size = 1.4))






# Figure 7: Correlation Analysis of Microbes and Metabolites from PLS-DA Integration


# IDENTIFIERS:
# Microbe identification table
tax <- asv_counts %>% 
        tax_table() %>% 
        data.frame() %>% 
        rownames_to_column('ASV')
# Metabolite identification table
mets <- featIDs %>% 
        dplyr::select(Compound, Accepted.Description, compound.class)

# Loading variables (from sPLS-DA)
micro.comp1 <- selectVar(diablo.model)$MIC$value %>% mutate(comp = 'comp1')
micro.comp2 <- selectVar(diablo.model, comp = 2)$MIC$value %>% mutate(comp = 'comp2')
metab.comp1 <- selectVar(diablo.model)$MET$value %>% mutate(comp = 'comp1')
metab.comp2 <- selectVar(diablo.model, comp = 2)$MET$value %>% mutate(comp = 'comp2')

mics <- rbind(micro.comp1, micro.comp2) %>% 
        rownames_to_column('ASV') %>% 
        left_join(tax) %>% 
        mutate(MInames = ifelse(is.na(Species), 
                                paste0(Genus, "_", ASV), 
                                paste0(Genus, "_", Species, "_", ASV))) 

metab <- rbind(metab.comp1, metab.comp2) %>% 
        rownames_to_column('Compound') %>% 
        left_join(mets) %>% 
        mutate(MBnames = paste0(Accepted.Description, ' ', '(', Compound, ")")) %>% 
        arrange(Compound)


# Spearman correlation
sub.multi <- multiset_data %>% 
        dplyr::select(SHIME, reactor, day, period, metab$Compound, mics$ASV) %>% 
        pivot_longer(starts_with('MB'), names_to = 'Compound', values_to = 'intensity') %>% 
        pivot_longer(starts_with('ASV'), names_to = 'ASV', values_to = 'clr') %>% 
        group_by(Compound, ASV) %>% 
        nest() %>% 
        mutate(sp.cor = purrr::map(data, function(x) cor.test(x$clr, x$intensity, method = 'spearman'))) %>% 
        mutate(rho = purrr::map_dbl(sp.cor, function(x) x$estimate[[1]]))

# Rho mat
spr.mat <- sub.multi %>% 
        ungroup() %>% 
        dplyr::select(-sp.cor, -data) %>% 
        left_join(mics %>% dplyr::select(ASV, MInames)) %>% 
        left_join(metab %>% dplyr::select(Compound, MBnames)) %>% 
        dplyr::select(-c(ASV, Compound)) %>% 
        pivot_wider(names_from = 'MInames', values_from = 'rho') %>% 
        column_to_rownames('MBnames')


# [save as 3000w x 1400h]
pheatmap::pheatmap(spr.mat, 
                   fontsize = 16,
                   angle_col = 315)







