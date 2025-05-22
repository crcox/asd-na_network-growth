library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
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

## no group ----
f0_nogroup <- update(fE, ~ . + RC1 + RC2 + RC3)
M0_nogroup <- netgrowr::mle_network_growth(
    f0_nogroup,
    data = modelvars_df,
    split_by = "vocab_step",
    label_with = "label"
)
saveRDS(M0_nogroup, "model-fits/baseline-model-nogroup.rds")


## no interaction ----
f0_nointeraction <- update(f0_nogroup, ~ . + group)
M0_nointeraction <- netgrowr::mle_network_growth(
    f0_nointeraction,
    data = modelvars_df,
    split_by = "vocab_step",
    label_with = "label"
)
saveRDS(M0_nointeraction, "model-fits/baseline-model-nointeraction.rds")


## with interaction ----
f0 <- update(f0_nogroup, ~ (.) * group)
M0 <- netgrowr::mle_network_growth(
    f0,
    data = modelvars_df,
    split_by = "vocab_step",
    label_with = "label"
)
saveRDS(M0, "model-fits/baseline-model.rds")


# Models with one growth value type ----
## no group ----
f1_nogroup <- list(
    assoc = list(
        pat = update(f0, ~ . + pat1_z),
        loa = update(f0, ~ . + loa1_z),
        pac = update(f0, ~ . + pac1_z)
    ),
    childes = list(
        pat = update(f0, ~ . + pat2_z),
        loa = update(f0, ~ . + loa2_z),
        pac = update(f0, ~ . + pac2_z)
    ),
    both = list(
        pat = update(f0, ~ . + pat1_z + pat2_z),
        loa = update(f0, ~ . + loa1_z + loa2_z),
        pac = update(f0, ~ . + pac1_z + pac2_z)
    )
)
M1_nogroup <- map_depth(f1_nogroup, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "step",
        label_with = "label"
    )
}, .progress = TRUE)
saveRDS(M1_nogroup, "model-fits/growth-models-nogroup.rds")


## no interaction ----
f1_nointeraction <- list(
    assoc = list(
        pat = update(f0, ~ . + pat1_z + group),
        loa = update(f0, ~ . + loa1_z + group),
        pac = update(f0, ~ . + pac1_z + group)
    ),
    childes = list(
        pat = update(f0, ~ . + pat2_z + group),
        loa = update(f0, ~ . + loa2_z + group),
        pac = update(f0, ~ . + pac2_z + group)
    ),
    both = list(
        pat = update(f0, ~ . + pat1_z + pat2_z + group ),
        loa = update(f0, ~ . + loa1_z + loa2_z + group ),
        pac = update(f0, ~ . + pac1_z + pac2_z + group )
    )
)
M1_nointeraction <- map_depth(f1_nointeraction, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "step",
        label_with = "label"
    )

}, .progress = TRUE)
saveRDS(M1_nointeraction, "model-fits/growth-models-nointeraction.rds")


## with interaction ----
f1 <- list(
    assoc = list(
        pat = update(f0, ~ (. + pat1_z) * group),
        loa = update(f0, ~ (. + loa1_z) * group),
        pac = update(f0, ~ (. + pac1_z) * group)
    ),
    childes = list(
        pat = update(f0, ~ (. + pat2_z) * group),
        loa = update(f0, ~ (. + loa2_z) * group),
        pac = update(f0, ~ (. + pac2_z) * group)
    ),
    both = list(
        pat = update(f0, ~ (. + pat1_z + pat2_z) * group ),
        loa = update(f0, ~ (. + loa1_z + loa2_z) * group ),
        pac = update(f0, ~ (. + pac1_z + pac2_z) * group )
    )
)
M1 <- map_depth(f1, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "step",
        label_with = "label"
    )
saveRDS(M1_nointeraction, "model-fits/growth-models.rds")

}, .progress = TRUE)
