library(dplyr)
library(purrr)
library(tidyr)
library(netgrowr)
library(parallel)
library(progressr)

# Load and prep data ----
modelvars_df <- readRDS("network/modelvars_vsoa_RC_z.rds")


# Empty (random selection) model ----
# +++ All words assigned equal probability
fE <- formula(vsoa_bin ~ 1)
ME <- netgrowr::mle_network_growth(fE, data = na.omit(modelvars_df), split_by = "vocab_step", label_with = "label")
saveRDS(ME, "model-fits/empty-model.rds")
# ME <- readRDS("model-fits/empty-model.rds")


# Psycholinguistic baseline models ----

# Principal Components Analysis
# Call: principal(r = select(mutate(confounds,
#                                   across(c(nphon, CHILDES_Freq),
#                                          list(log = ~log(.x + 1)))),
#                            nphon_log, CHILDES_Freq_log,
#                            BiphonProb.avg, PNDC.avg),
#                 nfactors = 3, rotate = "varimax")
# Standardized loadings (pattern matrix) based upon correlation matrix
#                    RC1   RC3  RC2   h2    u2 com
# nphon_log        -0.66 -0.56 0.21 0.79 0.206 2.2
# CHILDES_Freq_log  0.96  0.09 0.02 0.92 0.076 1.0
# BiphonProb.avg   -0.04 -0.03 0.99 0.99 0.011 1.0
# PNDC.avg          0.14  0.96 0.01 0.94 0.055 1.0
#
#                        RC1  RC3  RC2
# SS loadings           1.37 1.25 1.03
# Proportion Var        0.34 0.31 0.26
# Cumulative Var        0.34 0.65 0.91
# Proportion Explained  0.38 0.34 0.28
# Cumulative Proportion 0.38 0.72 1.00

## just group ----
f0_justgroup <- update(fE, ~ . + group)
M0_justgroup <- netgrowr::mle_network_growth(
    f0_justgroup,
    data = modelvars_df,
    split_by = "vocab_step",
    label_with = "label"
)
saveRDS(M0_justgroup, "model-fits/baseline-model-justgroup.rds")

## no group ----
f0_nogroup <- update(fE, ~ . + RC1 + RC2 + RC3)
M0_nogroup <- netgrowr::mle_network_growth(
    f0_nogroup,
    data = modelvars_df,
    split_by = "vocab_step",
    label_with = "label"
)
saveRDS(M0_nogroup, "model-fits/baseline-model-nogroup.rds")
# M0_nogroup <- readRDS("model-fits/baseline-model-nogroup.rds")


## no interaction ----
f0_nointeraction <- update(f0_nogroup, ~ . + group)
M0_nointeraction <- netgrowr::mle_network_growth(
    f0_nointeraction,
    data = modelvars_df,
    split_by = "vocab_step",
    label_with = "label"
)
saveRDS(M0_nointeraction, "model-fits/baseline-model-nointeraction.rds")
# M0_nointeraction <- readRDS("model-fits/baseline-model-nointeraction.rds")


## with interaction ----
f0 <- update(f0_nogroup, ~ (.) * group)
M0 <- netgrowr::mle_network_growth(
    f0,
    data = modelvars_df,
    split_by = "vocab_step",
    label_with = "label"
)
saveRDS(M0, "model-fits/baseline-model.rds")
M0 <- readRDS("model-fits/baseline-model.rds")



# Single RC models ----

## no group ----
fRC_1_nogroup <- list(
    RC1 = update(fE, ~ . + RC1),
    RC2 = update(fE, ~ . + RC2),
    RC3 = update(fE, ~ . + RC3)
)
MRC_1_nogroup <- map(fRC_1_nogroup,  ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "MRC 1 (no group)"))
saveRDS(MRC_1_nogroup, "model-fits/RC-1-model-nogroup.rds")


## no interaction ----
fRC_1_nointeraction <- list(
    RC1 = update(fE, ~ . + group + RC1),
    RC2 = update(fE, ~ . + group + RC2),
    RC3 = update(fE, ~ . + group + RC3)
)
MRC_1_nointeraction <- map(fRC_1_nointeraction,  ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "MRC 1 (no interaction)"))
saveRDS(MRC_1_nointeraction, "model-fits/RC-1-model-nointeraction.rds")


f0_nointeraction <- update(f0_nogroup, ~ . + group)
M0_nointeraction <- netgrowr::mle_network_growth(
    f0_nointeraction,
    data = modelvars_df,
    split_by = "vocab_step",
    label_with = "label"
)
saveRDS(M0_nointeraction, "model-fits/baseline-model-nointeraction.rds")
# M0_nointeraction <- readRDS("model-fits/baseline-model-nointeraction.rds")


## with interaction ----
fRC_1_full <- list(
    RC1 = update(fE, ~ (. + RC1) * group),
    RC2 = update(fE, ~ (. + RC2) * group),
    RC3 = update(fE, ~ (. + RC3) * group)
)
MRC_1_full <- map(fRC_1_full,  ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "MRC 1 (full)"))
saveRDS(MRC_1_full, "model-fits/RC-1-model-full.rds")


# Double RC models ----

## no group ----
fRC_2_nogroup <- list(
    RC1_RC2 = update(fE, ~ . + RC1 + RC2),
    RC1_RC3 = update(fE, ~ . + RC1 + RC3),
    RC2_RC3 = update(fE, ~ . + RC2 + RC3)
)
MRC_2_nogroup <- map(fRC_2_nogroup,  ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "MRC 2 (no group)"))
saveRDS(MRC_2_nogroup, "model-fits/RC-2-model-nogroup.rds")


## no interaction ----
fRC_2_nointeraction <- list(
    RC1_RC2 = update(fE, ~ . + group + RC1 + RC2),
    RC1_RC3 = update(fE, ~ . + group + RC1 + RC3),
    RC2_RC3 = update(fE, ~ . + group + RC2 + RC3)
)
MRC_2_nointeraction <- map(fRC_2_nointeraction,  ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "MRC 2 (no group)"))
saveRDS(MRC_2_nointeraction, "model-fits/RC-2-model-nointeraction.rds")

## with interaction ----
fRC_2_full <- list(
    RC1_RC2 = update(fE, ~ (. + RC1 + RC2) * group),
    RC1_RC3 = update(fE, ~ (. + RC1 + RC3) * group),
    RC2_RC3 = update(fE, ~ (. + RC2 + RC3) * group)
)
MRC_2_full <- map(fRC_2_full,  ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "MRC 2 (no group)"))
saveRDS(MRC_2_full, "model-fits/RC-2-model-full.rds")
