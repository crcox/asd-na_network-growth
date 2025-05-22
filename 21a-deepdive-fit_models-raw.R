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


# Models with one growth value type ----
## no group ----
f1_nogroup <- list(
    assoc = list(
        pat = update(f0_nogroup, ~ . + pat1_r),
        loa = update(f0_nogroup, ~ . + loa1_r),
        pac = update(f0_nogroup, ~ . + pac1_r)
    ),
    childes = list(
        pat = update(f0_nogroup, ~ . + pat2_r),
        loa = update(f0_nogroup, ~ . + loa2_r),
        pac = update(f0_nogroup, ~ . + pac2_r)
    ),
    both = list(
        pat = update(f0_nogroup, ~ . + pat1_r + pat2_r),
        loa = update(f0_nogroup, ~ . + loa1_r + loa2_r),
        pac = update(f0_nogroup, ~ . + pac1_r + pac2_r)
    )
)
M1_nogroup <- map_depth(f1_nogroup, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M1 (no group)"))
saveRDS(M1_nogroup, "model-fits/z-overall/growth-models-1-nogroup.rds")


## no interaction ----
f1_nointeraction <- list(
    assoc = list(
        pat = update(f0_nogroup, ~ . + pat1_r + group),
        loa = update(f0_nogroup, ~ . + loa1_r + group),
        pac = update(f0_nogroup, ~ . + pac1_r + group)
    ),
    childes = list(
        pat = update(f0_nogroup, ~ . + pat2_r + group),
        loa = update(f0_nogroup, ~ . + loa2_r + group),
        pac = update(f0_nogroup, ~ . + pac2_r + group)
    ),
    both = list(
        pat = update(f0_nogroup, ~ . + pat1_r + pat2_r + group ),
        loa = update(f0_nogroup, ~ . + loa1_r + loa2_r + group ),
        pac = update(f0_nogroup, ~ . + pac1_r + pac2_r + group )
    )
)
M1_nointeraction <- map_depth(f1_nointeraction, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M1 (no inter.)"))
saveRDS(M1_nointeraction, "model-fits/z-overall/growth-models-1-nointeraction.rds")


## with interaction ----
f1 <- list(
    assoc = list(
        pat = update(f0_nogroup, ~ (. + pat1_r) * group),
        loa = update(f0_nogroup, ~ (. + loa1_r) * group),
        pac = update(f0_nogroup, ~ (. + pac1_r) * group)
    ),
    childes = list(
        pat = update(f0_nogroup, ~ (. + pat2_r) * group),
        loa = update(f0_nogroup, ~ (. + loa2_r) * group),
        pac = update(f0_nogroup, ~ (. + pac2_r) * group)
    ),
    both = list(
        pat = update(f0_nogroup, ~ (. + pat1_r + pat2_r) * group ),
        loa = update(f0_nogroup, ~ (. + loa1_r + loa2_r) * group ),
        pac = update(f0_nogroup, ~ (. + pac1_r + pac2_r) * group )
    )
)
M1 <- map_depth(f1, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M1 (full)"))
saveRDS(M1, "model-fits/z-overall/growth-models-1-full.rds")


# Models with two growth value types ----
## no group ----
f2_nogroup <- list(
    assoc = list(
        pat = list(
            loa = update(f0_nogroup, ~ . + pat1_r + loa1_r),
            pac = update(f0_nogroup, ~ . + pat1_r + pac1_r)
        ),
        loa = list(
            pat = update(f0_nogroup, ~ . + loa1_r + pat1_r),
            pac = update(f0_nogroup, ~ . + loa1_r + pac1_r)
        ),
        pac = list(
            loa = update(f0_nogroup, ~ . + pat1_r + loa1_r),
            pac = update(f0_nogroup, ~ . + pat1_r + pac1_r)
        )
    ),
    childes = list(
        pat = list(
            loa = update(f0_nogroup, ~ . + pat2_r + loa2_r),
            pac = update(f0_nogroup, ~ . + pat2_r + pac2_r)
        ),
        loa = list(
            pat = update(f0_nogroup, ~ . + loa2_r + pat2_r),
            pac = update(f0_nogroup, ~ . + loa2_r + pac2_r)
        ),
        pac = list(
            loa = update(f0_nogroup, ~ . + pat2_r + loa2_r),
            pac = update(f0_nogroup, ~ . + pat2_r + pac2_r)
        )
    ),
    both = list(
        pat = list(
            loa = update(f0_nogroup, ~ . + pat1_r + pat2_r + loa1_r+ loa2_r),
            pac = update(f0_nogroup, ~ . + pat1_r + pat2_r + pac1_r+ pac2_r)
        ),
        loa = list(
            pat = update(f0_nogroup, ~ . + loa1_r + loa2_r + pat1_r + pat2_r),
            pac = update(f0_nogroup, ~ . + loa1_r + loa2_r + pac1_r + pac2_r)
        ),
        pac = list(
            loa = update(f0_nogroup, ~ . + pat1_r + pat2_r + loa1_r + loa2_r),
            pac = update(f0_nogroup, ~ . + pat1_r + pat2_r + pac1_r + pac2_r)
        )
    )
)
M2_nogroup <- map_depth(f2_nogroup, 3, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M2 (no group)"))
saveRDS(M2_nogroup, "model-fits/z-overall/growth-models-2-nogroup.rds")


## no interaction ----
f2_nointeraction <- list(
    assoc = list(
        pat = list(
            loa = update(f0_nogroup, ~ . + group + pat1_r + loa1_r),
            pac = update(f0_nogroup, ~ . + group + pat1_r + pac1_r)
        ),
        loa = list(
            pat = update(f0_nogroup, ~ . + group + loa1_r + pat1_r),
            pac = update(f0_nogroup, ~ . + group + loa1_r + pac1_r)
        ),
        pac = list(
            loa = update(f0_nogroup, ~ . + group + pat1_r + loa1_r),
            pac = update(f0_nogroup, ~ . + group + pat1_r + pac1_r)
        )
    ),
    childes = list(
        pat = list(
            loa = update(f0_nogroup, ~ . + group + pat2_r + loa2_r),
            pac = update(f0_nogroup, ~ . + group + pat2_r + pac2_r)
        ),
        loa = list(
            pat = update(f0_nogroup, ~ . + group + loa2_r + pat2_r),
            pac = update(f0_nogroup, ~ . + group + loa2_r + pac2_r)
        ),
        pac = list(
            loa = update(f0_nogroup, ~ . + group + pat2_r + loa2_r),
            pac = update(f0_nogroup, ~ . + group + pat2_r + pac2_r)
        )
    ),
    both = list(
        pat = list(
            loa = update(f0_nogroup, ~ . + group + pat1_r + pat2_r + loa1_r+ loa2_r),
            pac = update(f0_nogroup, ~ . + group + pat1_r + pat2_r + pac1_r+ pac2_r)
        ),
        loa = list(
            pat = update(f0_nogroup, ~ . + group + loa1_r + loa2_r + pat1_r + pat2_r),
            pac = update(f0_nogroup, ~ . + group + loa1_r + loa2_r + pac1_r + pac2_r)
        ),
        pac = list(
            loa = update(f0_nogroup, ~ . + group + pat1_r + pat2_r + loa1_r + loa2_r),
            pac = update(f0_nogroup, ~ . + group + pat1_r + pat2_r + pac1_r + pac2_r)
        )
    )
)
M2_nointeraction <- map_depth(f2_nointeraction, 3, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M2 (no interaction)"))
saveRDS(M2_nointeraction, "model-fits/z-overall/growth-models-2-nointeraction.rds")


## with interaction ----
f2 <- list(
    assoc = list(
        pat = list(
            loa = update(f0_nogroup, ~ (. + pat1_r + loa1_r) * group),
            pac = update(f0_nogroup, ~ (. + pat1_r + pac1_r) * group)
        ),
        loa = list(
            pat = update(f0_nogroup, ~ (. + loa1_r + pat1_r) * group),
            pac = update(f0_nogroup, ~ (. + loa1_r + pac1_r) * group)
        ),
        pac = list(
            loa = update(f0_nogroup, ~ (. + pat1_r + loa1_r) * group),
            pac = update(f0_nogroup, ~ (. + pat1_r + pac1_r) * group)
        )
    ),
    childes = list(
        pat = list(
            loa = update(f0_nogroup, ~ (. + pat2_r + loa2_r) * group),
            pac = update(f0_nogroup, ~ (. + pat2_r + pac2_r) * group)
        ),
        loa = list(
            pat = update(f0_nogroup, ~ (. + loa2_r + pat2_r) * group),
            pac = update(f0_nogroup, ~ (. + loa2_r + pac2_r) * group)
        ),
        pac = list(
            loa = update(f0_nogroup, ~ (. + pat2_r + loa2_r) * group),
            pac = update(f0_nogroup, ~ (. + pat2_r + pac2_r) * group)
        )
    ),
    both = list(
        pat = list(
            loa = update(f0_nogroup, ~ (. + pat1_r + pat2_r + loa1_r+ loa2_r) * group),
            pac = update(f0_nogroup, ~ (. + pat1_r + pat2_r + pac1_r+ pac2_r) * group)
        ),
        loa = list(
            pat = update(f0_nogroup, ~ (. + loa1_r + loa2_r + pat1_r + pat2_r) * group),
            pac = update(f0_nogroup, ~ (. + loa1_r + loa2_r + pac1_r + pac2_r) * group)
        ),
        pac = list(
            loa = update(f0_nogroup, ~ (. + pat1_r + pat2_r + loa1_r + loa2_r) * group),
            pac = update(f0_nogroup, ~ (. + pat1_r + pat2_r + pac1_r + pac2_r) * group)
        )
    )
)
M2 <- map_depth(f2, 3, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M2 (full)"))
saveRDS(M2_nogroup, "model-fits/z-overall/growth-models-2-full.rds")


# Models with three growth value types ----
## no group ----
f3_nogroup <- list(
    assoc = list(
        pat = list(
            loa = list(
                pac = update(f0_nogroup, ~ . + pat1_r + loa1_r + pac1_r)
            )
        )
    ),
    childes = list(
        pat = list(
            loa = list(
                pac = update(f0_nogroup, ~ . + pat2_r + loa2_r + pac2_r)
            )
        )
    ),
    both = list(
        pat = list(
            loa = list(
                pac = update(f0_nogroup, ~ . + pat1_r + loa1_r + pac1_r + pat2_r + loa2_r + pac2_r)
            )
        )
    )
)
M3_nogroup <- map_depth(f3_nogroup, 4, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M3 (no group)"))
saveRDS(M3_nogroup, "model-fits/z-overall/growth-models-3-nogroup.rds")


## no interaction ----
f3_nointeraction <- list(
    assoc = list(
        pat = list(
            loa = list(
                pac = update(f0_nogroup, ~ . + group + pat1_r + loa1_r + pac1_r)
            )
        )
    ),
    childes = list(
        pat = list(
            loa = list(
                pac = update(f0_nogroup, ~ . + group + pat2_r + loa2_r + pac2_r)
            )
        )
    ),
    both = list(
        pat = list(
            loa = list(
                pac = update(f0_nogroup, ~ . + group + pat1_r + loa1_r + pac1_r + pat2_r + loa2_r + pac2_r)
            )
        )
    )
)
M3_nointeraction <- map_depth(f3_nointeraction, 4, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M3 (no interaction)"))
saveRDS(M3_nointeraction, "model-fits/z-overall/growth-models-3-nointeraction.rds")


## with interaction ----
f3 <- list(
    assoc = list(
        pat = list(
            loa = list(
                pac = update(f0_nogroup, ~ (. + pat1_r + loa1_r + pac1_r) * group)
            )
        )
    ),
    childes = list(
        pat = list(
            loa = list(
                pac = update(f0_nogroup, ~ (. + pat2_r + loa2_r + pac2_r) * group)
            )
        )
    ),
    both = list(
        pat = list(
            loa = list(
                pac = update(f0_nogroup, ~ (. + pat1_r + loa1_r + pac1_r + pat2_r + loa2_r + pac2_r) * group)
            )
        )
    )
)
M3 <- map_depth(f3, 4, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M3 (full)"))
saveRDS(M3, "model-fits/z-overall/growth-models-3-full.rds")

