

# Analysis: Dynamic Regression modeling for n-of-1 Data
# Publication: Molecular Nutrition & Food Research Special Issue
# Project: SHIME
# Data: Short-chain fatty acids, Alpha diversity
# Author: Paige Jamieson
# Date: September 2024



# Environment ------------------------------------------------------------------

library(forecast)


# Short-chain fatty acids ----------------

# Read in SCFA Data
raw_fa <- read.csv("~/Documents/SHIME/Data/SHIME_SCFA_Results_mM.csv") %>% 
        column_to_rownames("sample.id") %>% 
        rename_with(~paste0("FA_", .x)) %>% 
        rownames_to_column("sample.id")

metadata <- read.csv("~/Documents/SHIME/Data/SHIME_metadata.csv")

totes <- left_join(raw_fa, metadata) %>% 
        column_to_rownames("sample.id") %>% 
        arrange(chron_order) %>% 
        filter(to_consider %in% c('Y', 'M', 'W')) %>% 
        mutate(shime_react = paste0(SHIME, "_", reactor)) %>% 
        mutate(reactor = factor(reactor, levels = c("AC", "TC", "DC"))) %>% 
        mutate(time = rep(seq(1, 17, 2), each = 6) %>% rep(2)) %>% 
        mutate(treatment = c(rep('A', each = 54), rep('B', each = 54)))

normtest <- fa_data %>% 
        group_by(SHIME, reactor, FA) %>% 
        nest() %>% 
        mutate(test = purrr::map(data, function(x) ad.test(x$conc))) %>% 
        mutate(pval = purrr::map_dbl(test, function(x) x$p.value))


nonnorm <- normtest %>% 
        dplyr::filter(pval < 0.05)
# Pass normality

# SCFA Data Long
fa_data <- totes %>% 
        pivot_longer(cols = c(FA_Acetic:FA_Butyric, FA_Total), names_to = 'FA', values_to = 'conc') %>% 
        dplyr::select(SHIME, reactor, time, treatment, FA, conc, day)



# Alpha Diversity ----------------------

# Read in phyloseq object
shime_ps <- readRDS("~/Documents/16S/SHIME/shime_ps.rds")

# Alpha Diversity

# Extract meta data
meta <- shime_ps %>% sample_data() %>% data.frame() %>% rownames_to_column('sample.id')

# Calculate Alpha diversity measures
adiv <- shime_ps %>%
        estimate_richness(measures = c('Observed', 'Shannon', 'Simpson')) %>% 
        rownames_to_column('sample.id') %>% 
        left_join(., meta) %>%
        dplyr::filter(period != 'stabilization') %>% 
        dplyr::filter(day != '31') %>% 
        mutate(time = rep(seq(1, 11, 2), each = 6) %>% rep(2)) %>% 
        mutate(dumby = c(rep(seq(1,9,2), each = 6), rep(seq(1,13,2), each = 6))) %>% 
        mutate(treatment = c(rep('A', each = 36), rep('B', each = 36)))

# Normality test
adivnorm <- adiv %>% 
        pivot_longer(cols = c(Observed, Shannon, Simpson), names_to = 'measure') %>% 
        group_by(measure, SHIME, period) %>% 
        nest() %>% 
        mutate(test = purrr::map(data, function(x) ad.test(x$value))) %>% 
        mutate(pval = purrr::map_dbl(test, function(x) x$p.value))
# Do not pass normality test

# ADIV Data Long
adiv_norm <- adiv %>% 
        mutate(across(c('Observed', 'Shannon', 'Simpson'), ~ log(.x))) %>% 
        pivot_longer(cols = c('Observed', 'Shannon', 'Simpson'), names_to = 'measure') %>% 
        dplyr::select(SHIME, reactor, time, treatment, measure, value, day, sampreact, dumby)







# SCFA Dynamic Regression Analysis ---------------------------------------------



#Visualize to determine time trends
ggplot(fa_data, aes(day, conc, group = shime_react, color = shime_react)) +
        geom_path() +
        facet_wrap(~FA, scales = 'free')


#Augmented Dickey-Fuller test (null hypothesis = series is non-stationary)
adf.data <- fa_data %>% 
        group_by(SHIME, reactor, FA) %>% 
        nest() %>% 
        mutate(test = purrr::map(data, function(x) adf.test(x$conc))) %>% 
        mutate(pval = purrr::map_dbl(test, function(x) round(x$p.value, 2)))
# All suggest non-stationarity





# ---------------------- ACETIC ACID ------------------------- #


stat.ace <- stat.data %>% 
        filter(FA == 'FA_Acetic') #%>% 
dplyr::select(SHIME, reactor, data, pval.trend, trend.mod, pval.mod)
# all have time trends

acetic <- fa_data %>% 
        filter(FA == 'FA_Acetic')


ggplot(acetic, aes(day, conc, group = shime_react, color = shime_react)) +
        geom_path() 


#Check autocorrelations


acetic %>% 
        dplyr::filter(SHIME == "B" & reactor == "DC") %>% 
        pull(conc) %>% 
        pacf(plot = T)


# SHIME A,AC: no autocorrelation
# SHIME A,TC, no autocorrelation
# SHIME A,DC, lag = 1
# SHIME B,AC, no autocorrelation
# SHIME B,TC, lag = 1, 4
# SHIME B,DC, lab = 1


# ------------ SHIME A - AC -------------#

acetic.aac <- acetic %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) %>% 
        mutate(L2 = dplyr::lag(conc, 2))

ace.aac.mod <- lm(conc ~ time + treatment + L1 + L2, acetic.aac)
summary(ace.aac.mod)
checkresiduals(ace.aac.mod)


ace.aac.coefs <- summary(ace.aac.mod)$coefficients


# ------------ SHIME A - TC -------------#

acetic.atc <- acetic %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'TC') 

ace.atc.mod <- lm(conc ~ time + treatment, acetic.atc)
summary(ace.atc.mod)
checkresiduals(ace.atc.mod)

ace.atc.coefs <- summary(ace.atc.mod)$coefficients

# ------------ SHIME A - DC -------------#

acetic.adc <- acetic %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(conc, 1))

ace.adc.mod <- lm(conc ~ time + treatment + L1, acetic.adc)
summary(ace.adc.mod)
checkresiduals(ace.adc.mod)

ace.adc.coefs <- summary(ace.adc.mod)$coefficients

# ------------ SHIME B - AC -------------#

acetic.bac <- acetic %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'AC')

ace.bac.mod <- lm(conc ~ time + treatment, acetic.bac)
summary(ace.bac.mod)
checkresiduals(ace.bac.mod)

ace.bac.coefs <- summary(ace.bac.mod)$coefficients

# ------------ SHIME B - TC -------------#

acetic.btc <- acetic %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) %>% 
        mutate(L4 = dplyr::lag(conc, 4))

ace.btc.mod <- lm(conc ~ time + treatment + L1 + L4, acetic.btc)
summary(ace.btc.mod)
checkresiduals(ace.btc.mod)

ace.btc.coefs <- summary(ace.btc.mod)$coefficients


# ------------ SHIME B - DC -------------#

acetic.bdc <- acetic %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(conc, 1))

ace.bdc.mod <- lm(conc ~ time + treatment + L1, acetic.bdc)
summary(ace.bdc.mod)
checkresiduals(ace.bac.mod)

ace.bdc.coefs <- summary(ace.bdc.mod)$coefficients


# Bind coefficients from all models
ace.coefs <- list(A.AC = ace.aac.coefs, 
                  A.TC = ace.atc.coefs,
                  A.DC = ace.adc.coefs,
                  B.AC = ace.bac.coefs,
                  B.TC = ace.btc.coefs,
                  B.DC = ace.bdc.coefs)


acetic.coefs <- do.call(rbind.data.frame, ace.coefs)


# ---------------------- BUTYRIC ACID ------------------------- #


stat.ace <- stat.data %>% 
        filter(FA == 'FA_Butyric') #%>% 
dplyr::select(SHIME, reactor, data, pval.trend, trend.mod, pval.mod)
# all have time trends

butyric <- fa_data %>% 
        filter(FA == 'FA_Butyric')



ggplot(butyric, aes(day, conc, group = shime_react, color = shime_react)) +
        geom_path() 


#Check autocorrelations

butyric %>% 
        dplyr::filter(SHIME == "B" & reactor == "DC") %>% 
        pull(conc) %>% 
        pacf(plot = T)


# SHIME A,AC: lag = 1
# SHIME A,TC, lag = 1
# SHIME A,DC, lag = 1
# SHIME B,AC, lag = 1
# SHIME B,TC, lag = 1, maybe 4
# SHIME B,DC, lag = 1


# ------------ SHIME A - AC -------------#

butyric.aac <- butyric %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

but.aac.mod <- lm(conc ~ time + treatment + L1, butyric.aac)
summary(but.aac.mod)
checkresiduals(but.aac.mod)

but.aac.coefs <- summary(but.aac.mod)$coefficients

# ------------ SHIME A - TC -------------#

butyric.atc <- butyric %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

but.atc.mod <- lm(conc ~ time + treatment + L1, butyric.atc)
summary(but.atc.mod)
checkresiduals(but.atc.mod)

but.atc.coefs <- summary(but.atc.mod)$coefficients

# ------------ SHIME A - DC -------------#

butyric.adc <- butyric %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

but.adc.mod <- lm(conc ~ time + treatment + L1, butyric.adc)
summary(but.adc.mod)
checkresiduals(but.adc.mod)

but.adc.coefs <- summary(but.adc.mod)$coefficients

# ------------ SHIME B - AC -------------#

butyric.bac <- butyric %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

but.bac.mod <- lm(conc ~ time + treatment + L1, butyric.bac)
summary(but.bac.mod)
checkresiduals(but.bac.mod)

but.bac.coefs <- summary(but.bac.mod)$coefficients


# ------------ SHIME B - TC -------------#

butyric.btc <- butyric %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) %>% 
        mutate(L4 = dplyr::lag(conc, 4))


but.btc.mod <- lm(conc ~ time + treatment + L1 + L4, butyric.btc)
summary(but.btc.mod)
checkresiduals(but.btc.mod)


but.btc.coefs <- summary(but.btc.mod)$coefficients


# ------------ SHIME B - DC -------------#

butyric.bdc <- butyric %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(conc, 1))

but.bdc.mod <- lm(conc ~ time + treatment + L1, butyric.bdc)
summary(but.bdc.mod)
checkresiduals(but.bdc.mod)

but.bdc.coefs <- summary(but.bdc.mod)$coefficients


# Bind coefficients from all models
but.coefs <- list(A.AC = but.aac.coefs, 
                  A.TC = but.atc.coefs,
                  A.DC = but.adc.coefs,
                  B.AC = but.bac.coefs,
                  B.TC = but.btc.coefs,
                  B.DC = but.bdc.coefs)


butyric.coefs <- do.call(rbind.data.frame, but.coefs)


# ---------------------- PROPIONIC ACID ------------------------- #


stat.ace <- stat.data %>% 
        filter(FA == 'FA_Propionic') #%>% 
dplyr::select(SHIME, reactor, data, pval.trend, trend.mod, pval.mod)
# all have time trends

propionic <- fa_data %>% 
        filter(FA == 'FA_Propionic')


ggplot(propionic, aes(day, conc, group = shime_react, color = shime_react)) +
        geom_path() 


#Check autocorrelations

propionic %>% 
        dplyr::filter(SHIME == "B" & reactor == "DC") %>% 
        pull(conc) %>% 
        pacf(plot = T)


# SHIME A,AC: lag = 1
# SHIME A,TC, lag = 1
# SHIME A,DC, lag = 1,2
# SHIME B,AC, no autocorrelation
# SHIME B,TC, 3,4
# SHIME B,DC, no autocorrelation


# ------------ SHIME A - AC -------------#

propionic.aac <- propionic %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

pro.aac.mod <- lm(conc ~ time + treatment + L1, propionic.aac)
summary(pro.aac.mod)
checkresiduals(pro.aac.mod)

pro.aac.coefs <- summary(pro.aac.mod)$coefficients


# ------------ SHIME A - TC -------------#

propionic.atc <- propionic %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

pro.atc.mod <- lm(conc ~ time + treatment + L1, propionic.atc)
summary(pro.atc.mod)
checkresiduals(pro.atc.mod)

pro.atc.coefs <- summary(pro.atc.mod)$coefficients

# ------------ SHIME A - DC -------------#

propionic.adc <- propionic %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

pro.adc.mod <- lm(conc ~ time + treatment + L1, propionic.adc)
summary(pro.adc.mod)
checkresiduals(pro.adc.mod)

pro.adc.coefs <- summary(pro.adc.mod)$coefficients


# ------------ SHIME B - AC -------------#

propionic.bac <- propionic %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

pro.bac.mod <- lm(conc ~ time + treatment, propionic.bac)
summary(pro.bac.mod)
checkresiduals(pro.bac.mod)

pro.bac.coefs <- summary(pro.bac.mod)$coefficients

# ------------ SHIME B - TC -------------#

propionic.btc <- propionic %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) %>% 
        mutate(L4 = dplyr::lag(conc, 4)) 

pro.btc.mod <- lm(conc ~ time + treatment + L4, propionic.btc)
summary(pro.btc.mod)
checkresiduals(pro.btc.mod)

pro.btc.coefs <- summary(pro.btc.mod)$coefficients

# ------------ SHIME B - DC -------------#

propionic.bdc <- propionic %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

pro.bdc.mod <- lm(conc ~ time + treatment, propionic.bdc)
summary(pro.bdc.mod)
checkresiduals(pro.bdc.mod)


pro.bdc.coefs <- summary(pro.bdc.mod)$coefficients


# Bind coefficients from all models
pro.coefs <- list(A.AC = pro.aac.coefs, 
                  A.TC = pro.atc.coefs,
                  A.DC = pro.adc.coefs,
                  B.AC = pro.bac.coefs,
                  B.TC = pro.btc.coefs,
                  B.DC = pro.bdc.coefs)


propionic.coefs <- do.call(rbind.data.frame, pro.coefs)




# ---------------------- TOTAL SCFA ------------------------- #


stat.ace <- stat.data %>% 
        filter(FA == "FA_Total") #%>% 
dplyr::select(SHIME, reactor, data, pval.trend, trend.mod, pval.mod)
# all have time trends

total <- fa_data %>% 
        filter(FA == 'FA_Total')


ggplot(total, aes(day, conc, group = shime_react, color = shime_react)) +
        geom_path() 


#Check autocorrelations

total %>% 
        dplyr::filter(SHIME == "B" & reactor == "DC") %>% 
        pull(conc) %>% 
        pacf(plot = T)


# SHIME A,AC: lag = 1
# SHIME A,TC, lag = 1
# SHIME A,DC, lag = 1
# SHIME B,AC, no autocorrelation
# SHIME B,TC, lag = 4
# SHIME B,DC, lag = 1


# ------------ SHIME A - AC -------------#

total.aac <- total %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

tot.aac.mod <- lm(conc ~ time + treatment + L1, total.aac)
summary(tot.aac.mod)
checkresiduals(tot.aac.mod)

tot.aac.coefs <- summary(tot.aac.mod)$coefficients

# ------------ SHIME A - TC -------------#

total.atc <- total %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

tot.atc.mod <- lm(conc ~ time + treatment + L1, total.atc)
summary(tot.atc.mod)
checkresiduals(tot.atc.mod)

tot.atc.coefs <- summary(tot.atc.mod)$coefficients

# ------------ SHIME A - DC -------------#

total.adc <- total %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

tot.adc.mod <- lm(conc ~ time + treatment + L1, total.adc)
summary(tot.adc.mod)
checkresiduals(tot.adc.mod)

tot.adc.coefs <- summary(tot.adc.mod)$coefficients

# ------------ SHIME B - AC -------------#

total.bac <- total %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

tot.bac.mod <- lm(conc ~ time + treatment, total.bac)
summary(tot.bac.mod)
checkresiduals(tot.bac.mod)


tot.bac.coefs <- summary(tot.bac.mod)$coefficients


# ------------ SHIME B - TC -------------#

total.btc <- total %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) %>% 
        mutate(L4 = dplyr::lag(conc, 4)) 

tot.btc.mod <- lm(conc ~ time + treatment + L4, total.btc)
summary(tot.btc.mod)
checkresiduals(tot.btc.mod)

tot.btc.coefs <- summary(tot.btc.mod)$coefficients


# ------------ SHIME B - DC -------------#

total.bdc <- total %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(conc, 1)) 

tot.bdc.mod <- lm(conc ~ time + treatment + L1, total.bdc)
summary(tot.bdc.mod)
checkresiduals(tot.bdc.mod)

tot.bdc.coefs <- summary(tot.bdc.mod)$coefficients



# Bind coefficients from all models
tot.coefs <- list(A.AC = tot.aac.coefs, 
                  A.TC = tot.atc.coefs,
                  A.DC = tot.adc.coefs,
                  B.AC = tot.bac.coefs,
                  B.TC = tot.btc.coefs,
                  B.DC = tot.bdc.coefs)


total.coefs <- do.call(rbind.data.frame, tot.coefs)





# Bind all results ---------------------------------

fa_results_list <- list(acetic = acetic.coefs,
                        propionic = propionic.coefs,
                        butyric = butyric.coefs,
                        total = total.coefs)

fa_results <- do.call(rbind.data.frame, fa_results_list)

#write.csv(fa_results, "~/Documents/SHIME/Data/fa_coefficients.csv")







# Alpha Diversity Dynammic Regression Analysis ---------------------------------


# ---------------------- Observed ------------------------- #

observed <- adiv_norm %>% 
        filter(measure == 'Observed')

ggplot(observed, aes(day, value, group = sampreact, color = sampreact)) +
        geom_path() 


#Check autocorrelations

observed %>% 
        dplyr::filter(SHIME == "B" & reactor == "DC") %>% 
        pull(value) %>% 
        pacf(plot = T)


# SHIME A,AC: no autocorrelation
# SHIME A,TC, lag = 1
# SHIME A,DC, no autocorrelation
# SHIME B,AC, no autocorrelation
# SHIME B,TC, no autocorrelation
# SHIME B,DC, no autocorrelation


# ------------ SHIME A - AC -------------#

observed.aac <- observed %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

obs.aac.mod <- lm(value ~ time + treatment, observed.aac)
summary(obs.aac.mod)
checkresiduals(obs.aac.mod)

obs.aac.coefs <- summary(obs.aac.mod)$coefficients

# ------------ SHIME A - TC -------------#

observed.atc <- observed %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

obs.atc.mod <- lm(value ~ time + treatment + L1, observed.atc)
summary(obs.atc.mod)
checkresiduals(obs.atc.mod)

obs.atc.coefs <- summary(obs.atc.mod)$coefficients

# ------------ SHIME A - DC -------------#

observed.adc <- observed %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

obs.adc.mod <- lm(value ~ time + treatment, observed.adc)
summary(obs.adc.mod)
checkresiduals(obs.adc.mod)

obs.adc.coefs <- summary(obs.adc.mod)$coefficients

# ------------ SHIME B - AC -------------#

observed.bac <- observed %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

obs.bac.mod <- lm(value ~ time + treatment, observed.bac)
summary(obs.bac.mod)
checkresiduals(obs.bac.mod)

obs.bac.coefs <- summary(obs.bac.mod)$coefficients


# ------------ SHIME B - TC -------------#

observed.btc <- observed %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

obs.btc.mod <- lm(value ~ time + treatment, observed.btc)
summary(obs.btc.mod)
checkresiduals(obs.btc.mod)

obs.btc.coefs <- summary(obs.btc.mod)$coefficients


# ------------ SHIME B - DC -------------#

observed.bdc <- observed %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

obs.bdc.mod <- lm(value ~ time + treatment, observed.bdc)
summary(obs.bdc.mod)
checkresiduals(obs.bdc.mod)

obs.bdc.coefs <- summary(obs.bdc.mod)$coefficients


# Bind coefficients from all models
obs.coefs <- list(A.AC = obs.aac.coefs, 
                  A.TC = obs.atc.coefs,
                  A.DC = obs.adc.coefs,
                  B.AC = obs.bac.coefs,
                  B.TC = obs.btc.coefs,
                  B.DC = obs.bdc.coefs)


observed.coefs <- do.call(rbind.data.frame, obs.coefs)






# ---------------------- Shannon ------------------------- #

shannon <- adiv_norm %>% 
        filter(measure == 'Shannon')

ggplot(shannon, aes(day, value, group = sampreact, color = sampreact)) +
        geom_path() 


#Check autocorrelations

shannon %>% 
        dplyr::filter(SHIME == "B" & reactor == "DC") %>% 
        pull(value) %>% 
        pacf(plot = T)


# SHIME A,AC: no autocorrelation
# SHIME A,TC, no autocorrelation
# SHIME A,DC, no autocorrelation
# SHIME B,AC, no autocorrelation
# SHIME B,TC, no autocorrelation
# SHIME B,DC, no autocorrelation


# ------------ SHIME A - AC -------------#

shannon.aac <- shannon %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

shan.aac.mod <- lm(value ~ time + treatment, shannon.aac)
summary(shan.aac.mod)
checkresiduals(shan.aac.mod)

shan.aac.coefs <- summary(shan.aac.mod)$coefficients

# ------------ SHIME A - TC -------------#

shannon.atc <- shannon %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

shan.atc.mod <- lm(value ~ time + treatment, shannon.atc)
summary(shan.atc.mod)
checkresiduals(shan.atc.mod)

shan.atc.coefs <- summary(shan.atc.mod)$coefficients

# ------------ SHIME A - DC -------------#

shannon.adc <- shannon %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

shan.adc.mod <- lm(value ~ time + treatment, shannon.adc)
summary(shan.adc.mod)
checkresiduals(shan.adc.mod)

shan.adc.coefs <- summary(shan.adc.mod)$coefficients

# ------------ SHIME B - AC -------------#

shannon.bac <- shannon %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

shan.bac.mod <- lm(value ~ time + treatment, shannon.bac)
summary(shan.bac.mod)
checkresiduals(shan.bac.mod)

shan.bac.coefs <- summary(shan.bac.mod)$coefficients

# ------------ SHIME B - TC -------------#

shannon.btc <- shannon %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

shan.btc.mod <- lm(value ~ time + treatment, shannon.btc)
summary(shan.btc.mod)
checkresiduals(shan.btc.mod)

shan.btc.coefs <- summary(shan.btc.mod)$coefficients

# ------------ SHIME B - DC -------------#

shannon.bdc <- shannon %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

shan.bdc.mod <- lm(value ~ time + treatment, shannon.bdc)
summary(shan.bdc.mod)
checkresiduals(shan.bdc.mod)

shan.bdc.coefs <- summary(shan.bdc.mod)$coefficients




# Bind coefficients from all models
shan.coefs <- list(A.AC = shan.aac.coefs, 
                   A.TC = shan.atc.coefs,
                   A.DC = shan.adc.coefs,
                   B.AC = shan.bac.coefs,
                   B.TC = shan.btc.coefs,
                   B.DC = shan.bdc.coefs)

shannon.coefs <- do.call(rbind.data.frame, shan.coefs)







# ---------------------- Simpson ------------------------- #

simpson <- adiv_norm %>% 
        filter(measure == 'Simpson')

ggplot(simpson, aes(day, value, group = sampreact, color = sampreact)) +
        geom_path() 


#Check autocorrelations

simpson %>% 
        dplyr::filter(SHIME == "B" & reactor == "DC") %>% 
        pull(value) %>% 
        pacf(plot = T)


# SHIME A,AC: no autocorrelation
# SHIME A,TC, no autocorrelation
# SHIME A,DC, no autocorrelation
# SHIME B,AC, no autocorrelation
# SHIME B,TC, no autocorrelation
# SHIME B,DC, no autocorrelation


# ------------ SHIME A - AC -------------#

simpson.aac <- simpson %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

simp.aac.mod <- lm(value ~ time + treatment, simpson.aac)
summary(simp.aac.mod)
checkresiduals(simp.aac.mod)

simp.aac.coefs <- summary(simp.aac.mod)$coefficients

# ------------ SHIME A - TC -------------#

simpson.atc <- simpson %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

simp.atc.mod <- lm(value ~ time + treatment, simpson.atc)
summary(simp.atc.mod)
checkresiduals(simp.atc.mod)

simp.atc.coefs <- summary(simp.atc.mod)$coefficients

# ------------ SHIME A - DC -------------#

simpson.adc <- simpson %>% 
        dplyr::filter(SHIME == 'A' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

simp.adc.mod <- lm(value ~ time + treatment, simpson.adc)
summary(simp.adc.mod)
checkresiduals(simp.adc.mod)

simp.adc.coefs <- summary(simp.adc.mod)$coefficients

# ------------ SHIME B - AC -------------#

simpson.bac <- simpson %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'AC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

simp.bac.mod <- lm(value ~ time + treatment, simpson.bac)
summary(simp.bac.mod)
checkresiduals(simp.bac.mod)

simp.bac.coefs <- summary(simp.bac.mod)$coefficients

# ------------ SHIME B - TC -------------#

simpson.btc <- simpson %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'TC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

simp.btc.mod <- lm(value ~ time + treatment, simpson.btc)
summary(simp.btc.mod)
checkresiduals(simp.btc.mod)

simp.btc.coefs <- summary(simp.btc.mod)$coefficients

# ------------ SHIME B - DC -------------#

simpson.bdc <- simpson %>% 
        dplyr::filter(SHIME == 'B' & reactor == 'DC') %>% 
        mutate(L1 = dplyr::lag(value, 1)) 

simp.bdc.mod <- lm(value ~ time + treatment, simpson.bdc)
summary(simp.bdc.mod)
checkresiduals(simp.bdc.mod)

simp.bdc.coefs <- summary(simp.bdc.mod)$coefficients



# Bind coefficients from all models
simp.coefs <- list(A.AC = simp.aac.coefs, 
                   A.TC = simp.atc.coefs,
                   A.DC = simp.adc.coefs,
                   B.AC = simp.bac.coefs,
                   B.TC = simp.btc.coefs,
                   B.DC = simp.bdc.coefs)

simpson.coefs <- do.call(rbind.data.frame, simp.coefs)





# Bind all results --------------------------------

adiv_results_list <- list(observed = observed.coefs,
                          shannon = shannon.coefs,
                          simpson = simpson.coefs)

adiv_results <- do.call(rbind.data.frame, adiv_results_list)


# write.csv(adiv_results, "~/Documents/SHIME/Data/adiv_coefficients.csv")







