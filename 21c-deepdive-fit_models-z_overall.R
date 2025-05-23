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
f0_full <- formula(vsoa_bin ~ (RC1 + RC2 + RC3) * group)


# Models with one growth value type ----
## no group ----
f1_nogroup <- list(
    assoc = list(
        pat = update(f0_nogroup, ~ . + pat1_Z),
        loa = update(f0_nogroup, ~ . + loa1_Z),
        pac = update(f0_nogroup, ~ . + pac1_Z)
    ),
    childes = list(
        pat = update(f0_nogroup, ~ . + pat2_Z),
        loa = update(f0_nogroup, ~ . + loa2_Z),
        pac = update(f0_nogroup, ~ . + pac2_Z)
    ),
    both = list(
        pat = update(f0_nogroup, ~ . + pat1_Z + pat2_Z),
        loa = update(f0_nogroup, ~ . + loa1_Z + loa2_Z),
        pac = update(f0_nogroup, ~ . + pac1_Z + pac2_Z)
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
        pat = update(f0_nogroup, ~ . + pat1_Z + group),
        loa = update(f0_nogroup, ~ . + loa1_Z + group),
        pac = update(f0_nogroup, ~ . + pac1_Z + group)
    ),
    childes = list(
        pat = update(f0_nogroup, ~ . + pat2_Z + group),
        loa = update(f0_nogroup, ~ . + loa2_Z + group),
        pac = update(f0_nogroup, ~ . + pac2_Z + group)
    ),
    both = list(
        pat = update(f0_nogroup, ~ . + pat1_Z + pat2_Z + group ),
        loa = update(f0_nogroup, ~ . + loa1_Z + loa2_Z + group ),
        pac = update(f0_nogroup, ~ . + pac1_Z + pac2_Z + group )
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


## Group interacts with baseline ----
f1_baseinteraction <- list(
    assoc = list(
        pat = update(f0_full, ~ . + pat1_Z),
        loa = update(f0_full, ~ . + loa1_Z),
        pac = update(f0_full, ~ . + pac1_Z)
    ),
    childes = list(
        pat = update(f0_full, ~ . + pat2_Z),
        loa = update(f0_full, ~ . + loa2_Z),
        pac = update(f0_full, ~ . + pac2_Z)
    ),
    both = list(
        pat = update(f0_full, ~ . + pat1_Z + pat2_Z),
        loa = update(f0_full, ~ . + loa1_Z + loa2_Z),
        pac = update(f0_full, ~ . + pac1_Z + pac2_Z)
    )
)
M1_baseinteraction <- map_depth(f1_baseinteraction, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M1 (base inter.)"))
saveRDS(M1_baseinteraction, "model-fits/z-overall/growth-models-1-baseinteraction.rds")


## with interaction ----
f1 <- list(
    assoc = list(
        pat = update(f0_nogroup, ~ (. + pat1_Z) * group),
        loa = update(f0_nogroup, ~ (. + loa1_Z) * group),
        pac = update(f0_nogroup, ~ (. + pac1_Z) * group)
    ),
    childes = list(
        pat = update(f0_nogroup, ~ (. + pat2_Z) * group),
        loa = update(f0_nogroup, ~ (. + loa2_Z) * group),
        pac = update(f0_nogroup, ~ (. + pac2_Z) * group)
    ),
    both = list(
        pat = update(f0_nogroup, ~ (. + pat1_Z + pat2_Z) * group ),
        loa = update(f0_nogroup, ~ (. + loa1_Z + loa2_Z) * group ),
        pac = update(f0_nogroup, ~ (. + pac1_Z + pac2_Z) * group )
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
        pat_loa = update(f0_nogroup, ~ . + pat1_Z + loa1_Z),
        pat_pac = update(f0_nogroup, ~ . + pat1_Z + pac1_Z),
        loa_pac = update(f0_nogroup, ~ . + loa1_Z + pac1_Z)
    ),
    childes = list(
        pat_loa = update(f0_nogroup, ~ . + pat2_Z + loa2_Z),
        pat_pac = update(f0_nogroup, ~ . + pat2_Z + pac2_Z),
        loa_pac = update(f0_nogroup, ~ . + loa2_Z + pac2_Z)
    ),
    both = list(
        pat_loa = update(f0_nogroup, ~ . + pat1_Z + pat2_Z + loa1_Z+ loa2_Z),
        pat_pac = update(f0_nogroup, ~ . + pat1_Z + pat2_Z + pac1_Z+ pac2_Z),
        loa_pac = update(f0_nogroup, ~ . + loa1_Z + loa2_Z + pac1_Z + pac2_Z)
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
saveRDS(M2_nogroup, "model-fits/z-overall/growth-models-2-nogroup.rds")


## no interaction ----
f2_nointeraction <- list(
    assoc = list(
        pat_loa = update(f0_nogroup, ~ . + group + pat1_Z + loa1_Z),
        pat_pac = update(f0_nogroup, ~ . + group + pat1_Z + pac1_Z),
        loa_pac = update(f0_nogroup, ~ . + group + loa1_Z + pac1_Z)
    ),
    childes = list(
        pat_loa = update(f0_nogroup, ~ . + group + pat2_Z + loa2_Z),
        pat_pac = update(f0_nogroup, ~ . + group + pat2_Z + pac2_Z),
        loa_pac = update(f0_nogroup, ~ . + group + loa2_Z + pac2_Z)
    ),
    both = list(
        pat_loa = update(f0_nogroup, ~ . + group + pat1_Z + pat2_Z + loa1_Z+ loa2_Z),
        pat_pac = update(f0_nogroup, ~ . + group + pat1_Z + pat2_Z + pac1_Z+ pac2_Z),
        loa_pac = update(f0_nogroup, ~ . + group + loa1_Z + loa2_Z + pac1_Z + pac2_Z)
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
saveRDS(M2_nointeraction, "model-fits/z-overall/growth-models-2-nointeraction.rds")


## Group interacts with baseline ----
f2_baseinteraction <- list(
    assoc = list(
        pat_loa = update(f0_full, ~ . + pat1_Z + loa1_Z),
        pat_pac = update(f0_full, ~ . + pat1_Z + pac1_Z),
        loa_pac = update(f0_full, ~ . + loa1_Z + pac1_Z)
    ),
    childes = list(
        pat_loa = update(f0_full, ~ . + pat2_Z + loa2_Z),
        pat_pac = update(f0_full, ~ . + pat2_Z + pac2_Z),
        loa_pac = update(f0_full, ~ . + loa2_Z + pac2_Z)
    ),
    both = list(
        pat_loa = update(f0_full, ~ . + pat1_Z + pat2_Z + loa1_Z + loa2_Z),
        pat_pac = update(f0_full, ~ . + pat1_Z + pat2_Z + pac1_Z + pac2_Z),
        loa_pac = update(f0_full, ~ . + loa1_Z + loa2_Z + pac1_Z + pac2_Z)
    )
)
M2_baseinteraction <- map_depth(f2_baseinteraction, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M2 (base interaction)"))
saveRDS(M2_baseinteraction, "model-fits/z-overall/growth-models-2-baseinteraction.rds")


## with interaction ----
f2 <- list(
    assoc = list(
        pat_loa = update(f0_nogroup, ~ (. + pat1_Z + loa1_Z) * group),
        pat_pac = update(f0_nogroup, ~ (. + pat1_Z + pac1_Z) * group),
        loa_pac = update(f0_nogroup, ~ (. + loa1_Z + pac1_Z) * group)
    ),
    childes = list(
        pat_loa = update(f0_nogroup, ~ (. + pat2_Z + loa2_Z) * group),
        pat_pac = update(f0_nogroup, ~ (. + pat2_Z + pac2_Z) * group),
        loa_pac = update(f0_nogroup, ~ (. + loa2_Z + pac2_Z) * group)
    ),
    both = list(
        pat_loa = update(f0_nogroup, ~ (. + pat1_Z + pat2_Z + loa1_Z + loa2_Z) * group),
        pat_pac = update(f0_nogroup, ~ (. + pat1_Z + pat2_Z + pac1_Z + pac2_Z) * group),
        loa_pac = update(f0_nogroup, ~ (. + loa1_Z + loa2_Z + pac1_Z + pac2_Z) * group)
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
saveRDS(M2, "model-fits/z-overall/growth-models-2-full.rds")


# Models with three growth value types ----
## no group ----
f3_nogroup <- list(
    assoc = list(
        pat_loa_pac = update(f0_nogroup, ~ . + pat1_Z + loa1_Z + pac1_Z)
    ),
    childes = list(
        pat_loa_pac = update(f0_nogroup, ~ . + pat2_Z + loa2_Z + pac2_Z)
    ),
    both = list(
        pat_loa_pac = update(f0_nogroup, ~ . + pat1_Z + loa1_Z + pac1_Z + pat2_Z + loa2_Z + pac2_Z)
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
saveRDS(M3_nogroup, "model-fits/z-overall/growth-models-3-nogroup.rds")


## no interaction ----
f3_nointeraction <- list(
    assoc = list(
        pat_loa_pac = update(f0_nogroup, ~ . + group + pat1_Z + loa1_Z + pac1_Z)
    ),
    childes = list(
        pat_loa_pac = update(f0_nogroup, ~ . + group + pat2_Z + loa2_Z + pac2_Z)
    ),
    both = list(
        pat_loa_pac = update(f0_nogroup, ~ . + group + pat1_Z + loa1_Z + pac1_Z + pat2_Z + loa2_Z + pac2_Z)
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
saveRDS(M3_nointeraction, "model-fits/z-overall/growth-models-3-nointeraction.rds")


## Group interacts with baseline ----
f3_baseinteraction <- list(
    assoc = list(
        pat_loa_pac = update(f0_full, ~ . + pat1_Z + loa1_Z + pac1_Z)
    ),
    childes = list(
        pat_loa_pac = update(f0_full, ~ . + pat2_Z + loa2_Z + pac2_Z)
    ),
    both = list(
        pat_loa_pac = update(f0_full, ~ . + pat1_Z + loa1_Z + pac1_Z + pat2_Z + loa2_Z + pac2_Z)
    )
)
M3_baseinteraction <- map_depth(f3_baseinteraction, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "vocab_step",
        label_with = "label"
    )
}, .progress = list(name = "M3 (no interaction)"))
saveRDS(M3_baseinteraction, "model-fits/z-overall/growth-models-3-baseinteraction.rds")


## with interaction ----
f3 <- list(
    assoc = list(
        pat_loa_pac = update(f0_nogroup, ~ (. + pat1_Z + loa1_Z + pac1_Z) * group)
    ),
    childes = list(
        pat_loa_pac = update(f0_nogroup, ~ (. + pat2_Z + loa2_Z + pac2_Z) * group)
    ),
    both = list(
        pat_loa_pac = update(f0_nogroup, ~ (. + pat1_Z + loa1_Z + pac1_Z + pat2_Z + loa2_Z + pac2_Z) * group)
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
saveRDS(M3, "model-fits/z-overall/growth-models-3-full.rds")

