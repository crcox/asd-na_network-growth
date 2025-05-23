library(dplyr)
library(purrr)
library(tidyr)
library(netgrowr)
library(parallel)
library(progressr)


# Load and prep data ----
modelvars_df <- readRDS("network/modelvars_vsoa_RC_z.rds")

## formula stem ----
f0_nogroup <- formula(vsoa_bin ~ RC1 + RC2 + RC3)


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
saveRDS(M1_nogroup, "model-fits/raw/growth-models-1-nogroup.rds")


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
saveRDS(M1_nointeraction, "model-fits/raw/growth-models-1-nointeraction.rds")


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
saveRDS(M1, "model-fits/raw/growth-models-1-full.rds")


# Models with two growth value types ----
## no group ----
f2_nogroup <- list(
    assoc = list(
        pat_loa = update(f0_nogroup, ~ . + pat1_r + loa1_r),
        pat_pac = update(f0_nogroup, ~ . + pat1_r + pac1_r),
        loa_pac = update(f0_nogroup, ~ . + loa1_r + pac1_r)
    ),
    childes = list(
        pat_loa = update(f0_nogroup, ~ . + pat2_r + loa2_r),
        pat_pac = update(f0_nogroup, ~ . + pat2_r + pac2_r),
        loa_pac = update(f0_nogroup, ~ . + loa2_r + pac2_r)
    ),
    both = list(
        pat_loa = update(f0_nogroup, ~ . + pat1_r + pat2_r + loa1_r + loa2_r),
        pat_pac = update(f0_nogroup, ~ . + pat1_r + pat2_r + pac1_r + pac2_r),
        loa_pac = update(f0_nogroup, ~ . + loa1_r + loa2_r + pac1_r + pac2_r)
    )
)
M2_nogroup <- map_depth(f2_nogroup, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M2 (no group)"))
saveRDS(M2_nogroup, "model-fits/raw/growth-models-2-nogroup.rds")


## no interaction ----
f2_nointeraction <- list(
    assoc = list(
        pat_loa = update(f0_nogroup, ~ . + group + pat1_r + loa1_r),
        pat_pac = update(f0_nogroup, ~ . + group + pat1_r + pac1_r),
        loa_pac = update(f0_nogroup, ~ . + group + loa1_r + pac1_r)
    ),
    childes = list(
        pat_loa = update(f0_nogroup, ~ . + group + pat2_r + loa2_r),
        pat_pac = update(f0_nogroup, ~ . + group + pat2_r + pac2_r),
        loa_pac = update(f0_nogroup, ~ . + group + loa2_r + pac2_r)
    ),
    both = list(
        pat_loa = update(f0_nogroup, ~ . + group + pat1_r + pat2_r + loa1_r+ loa2_r),
        pat_pac = update(f0_nogroup, ~ . + group + pat1_r + pat2_r + pac1_r+ pac2_r),
        loa_pac = update(f0_nogroup, ~ . + group + loa1_r + loa2_r + pac1_r + pac2_r)
    )
)
M2_nointeraction <- map_depth(f2_nointeraction, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M2 (no interaction)"))
saveRDS(M2_nointeraction, "model-fits/raw/growth-models-2-nointeraction.rds")


## with interaction ----
f2 <- list(
    assoc = list(
        pat_loa = update(f0_nogroup, ~ (. + pat1_r + loa1_r) * group),
        pat_pac = update(f0_nogroup, ~ (. + pat1_r + pac1_r) * group),
        loa_pac = update(f0_nogroup, ~ (. + loa1_r + pac1_r) * group)
    ),
    childes = list(
        pat_loa = update(f0_nogroup, ~ (. + pat2_r + loa2_r) * group),
        pat_pac = update(f0_nogroup, ~ (. + pat2_r + pac2_r) * group),
        loa_pac = update(f0_nogroup, ~ (. + loa2_r + pac2_r) * group)
    ),
    both = list(
        pat_loa = update(f0_nogroup, ~ (. + pat1_r + pat2_r + loa1_r+ loa2_r) * group),
        pat_pac = update(f0_nogroup, ~ (. + pat1_r + pat2_r + pac1_r+ pac2_r) * group),
        loa_pac = update(f0_nogroup, ~ (. + loa1_r + loa2_r + pac1_r + pac2_r) * group)
    )
)
M2 <- map_depth(f2, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M2 (full)"))
saveRDS(M2_nogroup, "model-fits/raw/growth-models-2-full.rds")


# Models with three growth value types ----
## no group ----
f3_nogroup <- list(
    assoc = list(
        pat_loa_pac = update(f0_nogroup, ~ . + pat1_r + loa1_r + pac1_r)
    ),
    childes = list(
        pat_loa_pac = update(f0_nogroup, ~ . + pat2_r + loa2_r + pac2_r)
    ),
    both = list(
        pat_loa_pac = update(f0_nogroup, ~ . + pat1_r + loa1_r + pac1_r + pat2_r + loa2_r + pac2_r)
    )
)
M3_nogroup <- map_depth(f3_nogroup, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M3 (no group)"))
saveRDS(M3_nogroup, "model-fits/raw/growth-models-3-nogroup.rds")


## no interaction ----
f3_nointeraction <- list(
    assoc = list(
        pat_loa_pac = update(f0_nogroup, ~ . + group + pat1_r + loa1_r + pac1_r)
    ),
    childes = list(
        pat_loa_pac = update(f0_nogroup, ~ . + group + pat2_r + loa2_r + pac2_r)
    ),
    both = list(
        pat_loa_pac = update(f0_nogroup, ~ . + group + pat1_r + loa1_r + pac1_r + pat2_r + loa2_r + pac2_r)
    )
)
M3_nointeraction <- map_depth(f3_nointeraction, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M3 (no interaction)"))
saveRDS(M3_nointeraction, "model-fits/raw/growth-models-3-nointeraction.rds")


## with interaction ----
f3 <- list(
    assoc = list(
        pat_loa_pac = update(f0_nogroup, ~ (. + pat1_r + loa1_r + pac1_r) * group)
    ),
    childes = list(
        pat_loa_pac = update(f0_nogroup, ~ (. + pat2_r + loa2_r + pac2_r) * group)
    ),
    both = list(
        pat_loa_pac = update(f0_nogroup, ~ (. + pat1_r + loa1_r + pac1_r + pat2_r + loa2_r + pac2_r) * group)
    )
)
M3 <- map_depth(f3, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M3 (full)"))
saveRDS(M3, "model-fits/raw/growth-models-3-full.rds")

