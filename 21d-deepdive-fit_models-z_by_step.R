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
        pat = update(f0_nogroup, ~ . + pat1_z),
        loa = update(f0_nogroup, ~ . + loa1_z),
        pac = update(f0_nogroup, ~ . + pac1_z)
    ),
    childes = list(
        pat = update(f0_nogroup, ~ . + pat2_z),
        loa = update(f0_nogroup, ~ . + loa2_z),
        pac = update(f0_nogroup, ~ . + pac2_z)
    ),
    both = list(
        pat = update(f0_nogroup, ~ . + pat1_z + pat2_z),
        loa = update(f0_nogroup, ~ . + loa1_z + loa2_z),
        pac = update(f0_nogroup, ~ . + pac1_z + pac2_z)
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
saveRDS(M1_nogroup, "model-fits/z-by-step/growth-models-1-nogroup.rds")


## no interaction ----
f1_nointeraction <- list(
    assoc = list(
        pat = update(f0_nogroup, ~ . + pat1_z + group),
        loa = update(f0_nogroup, ~ . + loa1_z + group),
        pac = update(f0_nogroup, ~ . + pac1_z + group)
    ),
    childes = list(
        pat = update(f0_nogroup, ~ . + pat2_z + group),
        loa = update(f0_nogroup, ~ . + loa2_z + group),
        pac = update(f0_nogroup, ~ . + pac2_z + group)
    ),
    both = list(
        pat = update(f0_nogroup, ~ . + pat1_z + pat2_z + group ),
        loa = update(f0_nogroup, ~ . + loa1_z + loa2_z + group ),
        pac = update(f0_nogroup, ~ . + pac1_z + pac2_z + group )
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
saveRDS(M1_nointeraction, "model-fits/z-by-step/growth-models-1-nointeraction.rds")


## Group interacts with baseline ----
f1_baseinteraction <- list(
    assoc = list(
        pat = update(f0_full, ~ . + pat1_z),
        loa = update(f0_full, ~ . + loa1_z),
        pac = update(f0_full, ~ . + pac1_z)
    ),
    childes = list(
        pat = update(f0_full, ~ . + pat2_z),
        loa = update(f0_full, ~ . + loa2_z),
        pac = update(f0_full, ~ . + pac2_z)
    ),
    both = list(
        pat = update(f0_full, ~ . + pat1_z + pat2_z),
        loa = update(f0_full, ~ . + loa1_z + loa2_z),
        pac = update(f0_full, ~ . + pac1_z + pac2_z)
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
saveRDS(M1_baseinteraction, "model-fits/z-by-step/growth-models-1-baseinteraction.rds")


## with interaction ----
f1 <- list(
    assoc = list(
        pat = update(f0_nogroup, ~ (. + pat1_z) * group),
        loa = update(f0_nogroup, ~ (. + loa1_z) * group),
        pac = update(f0_nogroup, ~ (. + pac1_z) * group)
    ),
    childes = list(
        pat = update(f0_nogroup, ~ (. + pat2_z) * group),
        loa = update(f0_nogroup, ~ (. + loa2_z) * group),
        pac = update(f0_nogroup, ~ (. + pac2_z) * group)
    ),
    both = list(
        pat = update(f0_nogroup, ~ (. + pat1_z + pat2_z) * group ),
        loa = update(f0_nogroup, ~ (. + loa1_z + loa2_z) * group ),
        pac = update(f0_nogroup, ~ (. + pac1_z + pac2_z) * group )
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
saveRDS(M1, "model-fits/z-by-step/growth-models-1-full.rds")


# Models with two growth value types ----
## no group ----
f2_nogroup <- list(
    assoc = list(
        pat_loa = update(f0_nogroup, ~ . + pat1_z + loa1_z),
        pat_pac = update(f0_nogroup, ~ . + pat1_z + pac1_z),
        loa_pac = update(f0_nogroup, ~ . + loa1_z + pac1_z)
    ),
    childes = list(
        pat_loa = update(f0_nogroup, ~ . + pat2_z + loa2_z),
        pat_pac = update(f0_nogroup, ~ . + pat2_z + pac2_z),
        loa_pac = update(f0_nogroup, ~ . + loa2_z + pac2_z)
    ),
    both = list(
        pat_loa = update(f0_nogroup, ~ . + pat1_z + pat2_z + loa1_z+ loa2_z),
        pat_pac = update(f0_nogroup, ~ . + pat1_z + pat2_z + pac1_z+ pac2_z),
        loa_pac = update(f0_nogroup, ~ . + loa1_z + loa2_z + pac1_z + pac2_z)
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
saveRDS(M2_nogroup, "model-fits/z-by-step/growth-models-2-nogroup.rds")


## no interaction ----
f2_nointeraction <- list(
    assoc = list(
        pat_loa = update(f0_nogroup, ~ . + group + pat1_z + loa1_z),
        pat_pac = update(f0_nogroup, ~ . + group + pat1_z + pac1_z),
        loa_pac = update(f0_nogroup, ~ . + group + loa1_z + pac1_z)
    ),
    childes = list(
        pat_loa = update(f0_nogroup, ~ . + group + pat2_z + loa2_z),
        pat_pac = update(f0_nogroup, ~ . + group + pat2_z + pac2_z),
        loa_pac = update(f0_nogroup, ~ . + group + loa2_z + pac2_z)
    ),
    both = list(
        pat_loa = update(f0_nogroup, ~ . + group + pat1_z + pat2_z + loa1_z + loa2_z),
        pat_pac = update(f0_nogroup, ~ . + group + pat1_z + pat2_z + pac1_z + pac2_z),
        loa_pac = update(f0_nogroup, ~ . + group + loa1_z + loa2_z + pac1_z + pac2_z)
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
saveRDS(M2_nointeraction, "model-fits/z-by-step/growth-models-2-nointeraction.rds")


## Group interacts with baseline ----
f2_baseinteraction <- list(
    assoc = list(
        pat_loa = update(f0_full, ~ . + pat1_z + loa1_z),
        pat_pac = update(f0_full, ~ . + pat1_z + pac1_z),
        loa_pac = update(f0_full, ~ . + loa1_z + pac1_z)
    ),
    childes = list(
        pat_loa = update(f0_full, ~ . + pat2_z + loa2_z),
        pat_pac = update(f0_full, ~ . + pat2_z + pac2_z),
        loa_pac = update(f0_full, ~ . + loa2_z + pac2_z)
    ),
    both = list(
        pat_loa = update(f0_full, ~ . + pat1_z + pat2_z + loa1_z + loa2_z),
        pat_pac = update(f0_full, ~ . + pat1_z + pat2_z + pac1_z + pac2_z),
        loa_pac = update(f0_full, ~ . + loa1_z + loa2_z + pac1_z + pac2_z)
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
saveRDS(M2_baseinteraction, "model-fits/z-by-step/growth-models-2-baseinteraction.rds")


## with interaction ----
f2 <- list(
    assoc = list(
        pat_loa = update(f0_nogroup, ~ (. + pat1_z + loa1_z) * group),
        pat_pac = update(f0_nogroup, ~ (. + pat1_z + pac1_z) * group),
        loa_pac = update(f0_nogroup, ~ (. + loa1_z + pac1_z) * group)
    ),
    childes = list(
        pat_loa = update(f0_nogroup, ~ (. + pat2_z + loa2_z) * group),
        pat_pac = update(f0_nogroup, ~ (. + pat2_z + pac2_z) * group),
        loa_pac = update(f0_nogroup, ~ (. + loa2_z + pac2_z) * group)
    ),
    both = list(
        pat_loa = update(f0_nogroup, ~ (. + pat1_z + pat2_z + loa1_z+ loa2_z) * group),
        pat_pac = update(f0_nogroup, ~ (. + pat1_z + pat2_z + pac1_z+ pac2_z) * group),
        loa_pac = update(f0_nogroup, ~ (. + loa1_z + loa2_z + pac1_z + pac2_z) * group)
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
saveRDS(M2_full, "model-fits/z-by-step/growth-models-2-full.rds")


# Models with three growth value types ----
## no group ----
f3_nogroup <- list(
    assoc = list(
        pat_loa_pac = update(f0_nogroup, ~ . + pat1_z + loa1_z + pac1_z)
    ),
    childes = list(
        pat_loa_pac = update(f0_nogroup, ~ . + pat2_z + loa2_z + pac2_z)
    ),
    both = list(
        pat_loa_pac = update(f0_nogroup, ~ . + pat1_z + loa1_z + pac1_z + pat2_z + loa2_z + pac2_z)
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
saveRDS(M3_nogroup, "model-fits/z-by-step/growth-models-3-nogroup.rds")


## no interaction ----
f3_nointeraction <- list(
    assoc = list(
        pat_loa_pac = update(f0_nogroup, ~ . + group + pat1_z + loa1_z + pac1_z)
    ),
    childes = list(
        pat_loa_pac = update(f0_nogroup, ~ . + group + pat2_z + loa2_z + pac2_z)
    ),
    both = list(
        pat_loa_pac = update(f0_nogroup, ~ . + group + pat1_z + loa1_z + pac1_z + pat2_z + loa2_z + pac2_z)
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
saveRDS(M3_nointeraction, "model-fits/z-by-step/growth-models-3-nointeraction.rds")


## Group interacts with baseline ----
f3_baseinteraction <- list(
    assoc = list(
        pat_loa_pac = update(f0_full, ~ . + pat1_z + loa1_z + pac1_z)
    ),
    childes = list(
        pat_loa_pac = update(f0_full, ~ . + pat2_z + loa2_z + pac2_z)
    ),
    both = list(
        pat_loa_pac = update(f0_full, ~ . + pat1_z + loa1_z + pac1_z + pat2_z + loa2_z + pac2_z)
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
saveRDS(M3_baseinteraction, "model-fits/z-by-step/growth-models-3-baseinteraction.rds")


## with interaction ----
f3 <- list(
    assoc = list(
        pat_loa_pac = update(f0_nogroup, ~ (. + pat1_z + loa1_z + pac1_z) * group)
    ),
    childes = list(
        pat_loa_pac = update(f0_nogroup, ~ (. + pat2_z + loa2_z + pac2_z) * group)
    ),
    both = list(
        pat_loa_pac = update(f0_nogroup, ~ (. + pat1_z + loa1_z + pac1_z + pat2_z + loa2_z + pac2_z) * group)
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
saveRDS(M3, "model-fits/z-by-step/growth-models-3-full.rds")

