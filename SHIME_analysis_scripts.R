


# Analysis scripts for publication in Molecular Nutrition & Food Research Special Issue
# Project: SHIME
# Data: XN Metabolites, Metabolomics, 16S rRNA sequencing
# Author: Paige Jamieson
# Date: September 2024


# Environment ------------------------------------------------------------------

library(phyloseq)
library(tidyverse)
library(nortest)
library(lmerTest)
library(lme4)
library(ggplot2)
library(microbiome)
library(mixOmics)
library(vegan)
library(pROC)
library(pheatmap)



# Functions to use throughout analyses:

# Funciton to calculate standard error
ser <- function(x){
        sd(x)/sqrt(length(x))
}

#Help functions:
log_helper <- function(x, min.val){
        log((x + sqrt(x^2 + min.val^2))/2)
}
#Pareto Scaling:
PS_helper <- function(x){
        (x - mean(x))/sqrt(sd(x, na.rm = T))
}	

#Transformation Functions:
#Log Scaling:
log_transform <- function(mtb){
        mtb_nz <- mtb[ ,which(apply(mtb, 2, sum) != 0)]
        min.val <- min(abs(mtb_nz[mtb_nz!=0]))/10
        mtb_log_trans <- apply(mtb_nz, 2, log_helper, min.val)
        return(mtb_log_trans)
}

#Pareto Scaling:
pareto_scale <- function(mtb){
        mtb_scaled <- apply(mtb, 2, PS_helper) 
        return(mtb_scaled)
}

#Function to find min value not zero
nzmin <- function(x){
        min(x[x > 0 ])
}

#Log2FC function to use on a dataframe
l2fcmin <- function(df, counts, var, condA, condB) {
        #Replace 0 values with 0.1% of min val
        df[[counts]] <- replace(df[[counts]], df[[counts]] <= 0, nzmin(df[[counts]]*0.001))
        c1 <- df[[counts]][sapply(df[[var]], function(x) x == condA)]
        c2 <- df[[counts]][sapply(df[[var]], function(x) x == condB)]
        l2fc <- mean(log2(c2[sapply(c2, is.finite)])) - mean(log2(c1[sapply(c1, is.finite)]))
        return(l2fc)
}



# Load Datasets ----------------------------------------------------

# Meta Data
meta<- read.csv("~/Documents/SHIME/Data/SHIME_metadata.csv") %>% 
        mutate(across(reactor, ~ factor(.x, levels = c('AC', 'TC', 'DC'))))

# XN Metabs
xn <- read.csv("~/Documents/SHIME/Data/SHIME_xn_metabs.csv") %>% 
        column_to_rownames("sample.id") %>% 
        rename_with(~paste0("XN_", .x)) %>% 
        rownames_to_column("sample.id")

# SCFA
scfa <- read.csv("~/Documents/SHIME/Data/SHIME_SCFA_Results_mM.csv") %>% 
        column_to_rownames("sample.id") %>% 
        rename_with(~paste0("FA_", .x)) %>% 
        rownames_to_column("sample.id")

# Load phyloseq object
shime_ps <- readRDS("~/Documents/16S/SHIME/shime_ps.rds") # ASV # = 718



# XN Metabolite Analysis ---------------------------------------------

# Join data
xn_combined <- right_join(meta, xn)

# Unpaired ttest considered by week, pvalues adjusted via BH
xn_test <- xn_combined %>% 
        filter(!week %in% c(3, 4)) %>% 
        pivot_longer(cols = c(starts_with('XN'), -XN_XN), names_to = 'metabolite', values_to = 'conc') %>% 
        group_by(metabolite, reactor, week) %>% 
        nest() %>% 
        mutate(ttest = purrr::map(data, function(x) t.test(x$conc ~ x$SHIME))) %>% 
        mutate(pval = purrr::map_dbl(ttest, function(x) x$p.value)) %>% 
        ungroup(reactor) %>% 
        mutate(padj = p.adjust(pval, method = 'BH'))




# Microbiome Analysis -----------------------------------------------

# Load phyloseq object
shime_ps <- readRDS("~/Documents/16S/SHIME/shime_ps.rds") # ASV # = 718  

# Create placeholder names for un-annotated genera using family
renames <- rownames(tax_table(shime_ps)[is.na(tax_table(shime_ps)[, 'Genus'])])
taxdf <- tax_table(shime_ps)[renames,]
renamed_genus <- unname(sapply(taxa_names(taxdf), function(x) paste0('f_', taxdf[x, 'Family'])))
# Add Family placehold names
tax_table(shime_ps)[renames, 'Genus'] <- renamed_genus
# Create placeholder names for un-annotated genera using order
rename_order <- tax_table(shime_ps) %>% data.frame() %>% filter(., grepl("f_NA", Genus)) %>% rownames()
taxadf <- tax_table(shime_ps)[rename_order,]
renamed_from_order <- unname(sapply(taxa_names(taxadf), function(x) paste0('o_', taxadf[x, 'Order'])))
# Add Order placeholder names 
tax_table(shime_ps)[rename_order, 'Genus'] <- renamed_from_order


# Filter out taxa seen fewer than 10 times in less than 10% of samples
asv_counts <- shime_ps %>% filter_taxa(function(x) sum(x > 10) > (0.1*length(x)), TRUE)
#Filter out low abundance (>1e-5) taxa >>> does not change number of ASVs 
counts_filtered <- asv_counts %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE) 
# Transform to relative abundance
asv_relab <- asv_counts %>% transform_sample_counts(function(x) x / sum(x))




# Beta Diversity Analysis (Aitchison Distance) ----------------------

# Subset data for control & treatment weeks to compare
counts_filtered <- asv_counts %>% 
        subset_samples(week != '5' & !period %in% c('stabilization', 'wash out'))

# CLR transform
clr <- counts_filtered %>% 
        rarefy_even_depth(rngseed = 44) %>% 
        otu_table() %>% 
        as.data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1)

# Extract meta data
meta <- counts_filtered %>% 
        sample_data() %>% 
        data.frame() %>% 
        rownames_to_column('sample.id') %>% 
        dplyr::filter(!period == 'stabilization') %>% 
        mutate(treatment = ifelse(period == 'control', 'pre XN', 'post XN'))

# Calculate Euclidean distances
euclid_dist <- clr %>% 
        stats::prcomp() 

# Extracts components
comps <- euclid_dist$x %>% 
        as.data.frame %>% 
        rownames_to_column('sample.id') %>% 
        right_join(meta)

# Plot PCA
aitch_pca <- ggplot(comps, aes(x = PC1, y = PC2, color = period, shape = SHIME))+
        geom_point(aes(size=3))




# PERMANOVA and Distance Based Redundancy Analysis (dbRDA) ----------

# CLR Transform OTU table in Phyloseq object
ps_clr <- microbiome::transform(counts_filtered, transform = 'clr')
# Calculate ordination
ordRDA <- ordinate(ps_clr, method = 'RDA', distance = 'euclidean')
# Plot
plot_ordination(ps_clr, ordRDA, color = 'period') + 
        geom_point(size = 1.5)


# Create distance matrix
distmat <- phyloseq::distance(ps_clr, method = 'euclidean')
# Extract metadata
mdata <- data.frame(sample_data(ps_clr))
# Verify variance is equal between groups
bdisp <- betadisper(distmat, mdata$period)
# Evaluate using a permutation test
permutest(bdisp) # p = 0.004 (period) (meaning any significance detected by 
#PERMANOVA may mean observed differences between groups are driven by spread rather than treatment effect)
adonis2(distmat ~ period, data = mdata)
# p-value = 0.008

# RDA analysis by treatment period
ordRDA <- ordinate(ps_clr, method = 'RDA', distance = 'euclidean', ~ period)

# Plot constrained analysis
plot_ordination(ps_clr, ordRDA, color = 'period') +
        geom_point(size = 1.5)
# treatment explains about 9.4% of the variance



# Differential Abundance Analysis ------------------------------------

cts_rare <- counts_filtered %>% 
        rarefy_even_depth(rngseed = 23) %>% 
        psmelt()

# Subset taxonomy data
tax <- tax_table(counts_filtered) %>% data.frame() %>% rownames_to_column('OTU') 

# Negative Binomial GLMM
glmDA <- cts_rare %>% 
        group_by(SHIME, OTU) %>% 
        nest() %>% 
        mutate(zpro = purrr::map_dbl(data, function(x) (sum(x$Abundance == 0) / length(x$Abundance))  )) %>%
        mutate(nbmod = purrr::map(data, function(x) try(glmer.nb(Abundance ~ period + (1|reactor), data = x, na.action = na.exclude)))) %>% 
        mutate(err = purrr::map_lgl(nbmod, function(x) any(class(x) == 'try-error'))) %>% 
        # Filter ASVs that produced model error
        filter(err == FALSE) %>% 
        # Filter ASVs with proportion of zeros > 85%
        dplyr::filter(zpro < 0.85) %>% 
        mutate(summ = purrr::map(nbmod, function(x) summary(x))) %>% 
        mutate(pval = purrr::map_dbl(summ, function(x) x$coefficients[[2,4]])) %>% 
        ungroup() %>% 
        mutate(padj = p.adjust(pval, method = 'BH')) 

# Filter for significance
glm.fil <- glmDA %>% 
        dplyr::select(OTU, SHIME, summ, padj) %>% 
        mutate(trt = purrr::map_dbl(summ, function(x) x$coefficients[[2,1]])) %>% 
        filter(padj <= 0.05 & padj > 0) %>% 
        left_join(tax, by = 'OTU')


# Determine number of taxa abundance that change with treatment
glm.fil %>% filter(SHIME == 'A') %>% filter(trt > 0 ) %>% dim() # 20 ASVs inc with trt
glm.fil %>% filter(SHIME == 'A') %>% filter(trt < 0 ) %>% dim() # 16 ASVs dec with trt
glm.fil %>% filter(SHIME == 'B') %>% filter(trt > 0 ) %>% dim() # 10 ASVs inc with trt
glm.fil %>% filter(SHIME == 'B') %>% filter(trt < 0 ) %>% dim() # 12 ASVs dec with trt





# Metabolomics ----------------------------------------------------



# Read in Progenesis QI datasets

# Positive ion mode
pos <- read.csv("~/Documents/SHIME/Metabolomics/SHIME_POS_092424.csv", skip = 2, na.string = c(""))
# Negative ion mode
neg <- read.csv("~/Documents/SHIME/Metabolomics/SHIME_NEG_092424.csv", skip = 2, na.strings = c(""))

# Read in sample identifiers
posids <- read.csv("~/Documents/SHIME/Data/pos_sample_identifiers.csv")
negids <- posids %>% mutate(run.id = gsub("pos", 'neg', run.id))

# Read in meta data
meta <- read.csv("~/Documents/SHIME/Data/SHIME_metadata.csv") %>% 
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
posQC <- read.csv("~/Documents/SHIME/Metabolomics/SHIME_POS_QC_092424.csv", skip = 2, na.string = c(""))
negQC <- read.csv("~/Documents/SHIME/Metabolomics/SHIME_NEG_QC_092424.csv", skip = 2, na.string = c(""))



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


ttest_helper <- function(dat, condition, measure, condA = 'control', condB = 'treatment'){
        cA <- dat[[measure]][sapply(dat[[condition]], function(x) x == condA)]
        cB <- dat[[measure]][sapply(dat[[condition]], function(x) x == condB)]
        M <- data.frame(ctrl = cA, trt = cB)
        t.test(Pair(trt, ctrl) ~ 1, data = M)
        
}


# Pairwise ttest before and after XN treatment
pairwise_by_reactor <- total_scaled %>% 
        column_to_rownames('sample.id') %>% 
        dplyr::select(-date, -day, -chron_order, -to_consider, -week) %>% 
        pivot_longer(starts_with('MB'), names_to = 'Compound', values_to = 'intensity') %>% 
        group_by(SHIME, reactor, Compound) %>% 
        nest() %>%
        mutate(ttest = purrr::map(data, function(x) ttest_helper(dat=x, 'period', 'intensity'))) %>% 
        mutate(pval = purrr::map_dbl(ttest, function(x) x$p.value)) %>% 
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
        mutate(across(A_AC:B_DC, ~replace_na(., 0 ))) %>%
        mutate(across(A_AC:B_DC, ~ p.adjust(.x, method = 'BH')))



# Multi-omic Integration Analysis (PLS-DA) -----------------------


# Prep Data for Multi-block PLS-DA 

# Pull metabolomics dataset
metab_scaled <- total_scaled %>% 
        dplyr::select(-c(date, chron_order, to_consider)) %>% 
        dplyr::select(-'MB_22.69_263.0925m/z', -'MB_21.08_173.1179m/z')

# Pull microbiome dataset
counts_filtered <- asv_counts %>% 
        subset_samples(week != '5' & !period %in% c('stabilization', 'wash out'))


# CLR transform
asvclr <- counts_filtered %>% 
        rarefy_even_depth(rngseed = 44) %>% 
        otu_table() %>% 
        data.frame() %>% 
        decostand(method = 'clr', pseudocount = 1) %>% 
        rownames_to_column('sample.id')
asvmd <- counts_filtered %>% 
        sample_data() %>% 
        data.frame() %>%
        rownames_to_column('sample.id') %>% 
        dplyr::select(-c(date, chron_order, to_consider, sampling.ID, samptreat, sampreact, phase))

# Join meta data
micdata <- left_join(asvmd, asvclr)


# Join microbiome and metabolome datasets so all rows cleanly match
multiset_data <- left_join(metab_scaled, micdata) %>% 
        column_to_rownames('sample.id') %>% 
        mutate(ind_fct = paste0(SHIME, '_', reactor)) %>% 
        mutate(DA = paste0(SHIME, "_", "period"))



# Initial sPLS analyses to determine correlation structure

# set a list of all the X dataframes
Xlist = list(MIC = multiset_data %>% dplyr::select(starts_with("ASV")),
             MET = multiset_data %>% dplyr::select(starts_with("MB")))
lapply(Xlist, dim)

# set treatment as the Y variable
trt <- as.factor(multiset_data$period)


list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(25, 25)


# MICRO and METAB
pls1 <- spls(Xlist[["MIC"]], Xlist[["MET"]], 
             keepX = list.keepX, keepY = list.keepY, 
             multilevel = multiset_data$ind_fct) 

# plot features of first PLS
plotVar(pls1, cutoff = 0.5, title = "(a) micro vs metabolome", 
        legend = c("micro", "metab"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# calculate correlation of micro and metabolome
cor(pls1$variates$X, pls1$variates$Y) # 0.972



# Multi-block PLS-DA Analysis ---------------------------


Xlist = list(MIC = multiset_data %>% dplyr::select(starts_with("ASV")),
             MET = multiset_data %>% dplyr::select(starts_with("MB")))

trt <- as.factor(multiset_data$period)


#Create the design matrix:
design <- matrix(0.5, nrow = 3, ncol = 3)
diag(design) <- 0


#Adjust our data to handle repeated measures
Cov <- multiset_data[, 'ind_fct', drop = F]
Xnew <- lapply(Xlist, function(i) withinVariation(X = i, design = Cov))



keepX <- list(MIC = rep(25, 2), MET = rep(25, 2))
# DIABLO model - with multilevel decomposition
diablo.model = block.splsda(X = Xnew, Y = trt,, keepX = keepX, ncomp = 2, design = design) 



#Create the consensus space by finding the average loadings of each compononet
conDiablo <- data.frame(Reduce('+', diablo.model$variates)/2) %>%
        cbind(diablo.model$Y) %>% 
        rename('condition' = 'diablo.model$Y') %>% 
        modify_at('treatment', as.character)


# Palette for plot
pal <- RColorBrewer::brewer.pal(4, 'Dark2')
names(pal) <- c('control', 'treatment')

# Plot consensus plot
conPlot <- ggplot(conDiablo, aes(comp1, comp2, color = condition)) +
        geom_point(size = 2) + 
        stat_ellipse(aes(group = condition, fill = condition), geom = 'polygon', alpha = 0.2) +
        theme(legend.position = 'bottom', legend.box = 'horizontal', 
              legend.title = element_blank()) +
        cowplot::theme_cowplot() +
        labs(x = 'Component 1', y = 'Component 2') +
        scale_color_manual(values = pal, labels = c('control', 'treatment'),
                           name = 'Condition') + 
        scale_fill_manual(values = pal, labels = c('control', 'treatment'),
                          name = 'Condition')


plotLoadings(diablo.model, comp = 1)
plotLoadings(diablo.model, comp = 2)



# Leave one out cross validation
perf.diablo = perf(diablo.model, validation = 'loo', auc = T, progressBar = F) 


# Extract out the data to make ROC curves
mic1 <- as.data.frame(perf.diablo$predict$nrep1$MIC$comp1)
mic2 <- as.data.frame(perf.diablo$predict$nrep1$MIC$comp2)
met1 <- as.data.frame(perf.diablo$predict$nrep1$MET$comp1)
met2 <- as.data.frame(perf.diablo$predict$nrep1$MET$comp2)


#Make the class matrix
classmat <- data.frame(class = diablo.model$Y) %>%
        mutate(ctrl = ifelse(class == 'control', 1, 0)) %>%
        mutate(trt = ifelse(class == 'treatment', 1, 0)) %>%
        dplyr::select(-class)


#Build a list of the outputs 
outputs <- list(micro_comp1 = mic1, micro_comp2 = mic2, metab_comp1 = met1, metab_comp2 = met2)


#Custom function to pull the classes and see how they do
pull_roc <- function(prediction_list, classes){
        purrr::map(prediction_list, function(x){
                map2(classes, x, roc)
        })
}




#Pull out the ROCS
rocs <- pull_roc(outputs, classmat)
#Make them into nice plots using ggroc
rocplots <- purrr::map(rocs, ggroc)
#Extract out the AUCs 
aucs <- purrr::map(rocs, function(x) purrr::map(x, auc))
pal <- RColorBrewer::brewer.pal(4, 'Dark2')
pal <- c('#E33255', '#3BE1AF')
names(pal) <- c('control', 'treatment')
palmb <- pal
names(palmb) <- c('ctrl', 'trt')



#Make ROC plots 
p1 <- rocplots[[1]] +
        geom_path(size = 2) + 
        geom_abline(intercept = 1) +
        cowplot::theme_cowplot() +
        ggtitle('Microbiome - Component 1')  +
        theme(legend.position = c(0.5,0.25)) +
        scale_color_manual(name = 'Condition',
                           values = palmb)

p2 <- rocplots[[2]] +
        geom_path(size = 2) + 
        geom_abline(intercept = 1) +
        cowplot::theme_cowplot() +
        ggtitle('Microbiome - Component 2')  +
        theme(legend.position = c(0.5,0.25)) +
        scale_color_manual(name = 'Condition',
                           values = palmb)

p3 <- rocplots[[3]] +
        geom_path(size = 2) + 
        geom_abline(intercept = 1) +
        cowplot::theme_cowplot() +
        ggtitle('Metabolome - Component 1')  +
        theme(legend.position = c(0.5,0.25)) +
        scale_color_manual(name = 'Condition',
                           values = palmb)

p4 <- rocplots[[4]] +
        geom_path(size = 2) + 
        geom_abline(intercept = 1) +
        cowplot::theme_cowplot() +
        ggtitle('Metabolome - Component 2')  +
        theme(legend.position = c(0.5,0.25)) +
        scale_color_manual(name = 'Condition',
                           values = palmb)

allrocs <- plot_grid(p1, p2, p3, p4, ncol = 2)





# Spearman correlation with final variables -----------------------------------

# Microbe identification table
tax <- asv_counts %>% 
        tax_table() %>% 
        data.frame() %>% 
        rownames_to_column('ASV')

# Metabolite identification table
mets <- featIDs %>% 
        dplyr::select(Compound, Accepted.Description, compound.class)

# Loading variables (from sPLS-DA)
micro.comp1 <- selectVar(diablo.model)$MIC$value %>% mutate(comp = 'comp1') %>% rownames_to_column('ASV')

micro.comp2 <- selectVar(diablo.model, comp = 2)$MIC$value %>% mutate(comp = 'comp2')%>% rownames_to_column('ASV')


subs <- micro.comp1$ASV %in% micro.comp2$ASV

# ASV repeats on compp2
to_rm <- micro.comp1$ASV[subs]


micro.comp2 <- selectVar(diablo.model, comp = 2)$MIC$value %>% mutate(comp = 'comp2')%>% rownames_to_column('ASV') %>% 
        dplyr::filter(!ASV %in% to_rm)
metab.comp1 <- selectVar(diablo.model)$MET$value %>% mutate(comp = 'comp1')
metab.comp2 <- selectVar(diablo.model, comp = 2)$MET$value %>% mutate(comp = 'comp2')


mics <- rbind(micro.comp1, micro.comp2) %>% 
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
        dplyr::select(SHIME, reactor, day, period, metab$Compound, unique(mics$ASV)) %>% 
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
        left_join(mics %>% dplyr::select(ASV, MInames), by = 'ASV') %>% 
        left_join(metab %>% dplyr::select(Compound, MBnames)) %>% 
        dplyr::select(-c(ASV, Compound)) %>% 
        pivot_wider(names_from = 'MInames', values_from = 'rho') %>% 
        column_to_rownames('MBnames')


# [save as 3000w x 1400h]
pheatmap::pheatmap(spr.mat, 
                   fontsize = 20,
                   angle_col = 315, 
                   border_color = 'grey90',
                   cellwidth = 30,
                   cellheight = 30)




