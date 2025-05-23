library(dplyr)
library(purrr)
library(tidyr)
library(netgrowr)


M1_nogroup_r <- readRDS("model-fits/raw/growth-models-1-nogroup-formula.rds")
M1_nointeraction_r <- readRDS("model-fits/raw/growth-models-1-nointeraction-formula.rds")
M1_full_r <- readRDS("model-fits/raw/growth-models-1-full-formula.rds")

M1_nogroup_Z <- readRDS("model-fits/z-overall/growth-models-1-nogroup-formula.rds")
M1_nointeraction_Z <- readRDS("model-fits/z-overall/growth-models-1-nointeraction-formula.rds")
M1_full_Z <- readRDS("model-fits/z-overall/growth-models-1-full-formula.rds")

M1_nogroup_z <- readRDS("model-fits/z-by-step/growth-models-1-nogroup-formula.rds")
M1_nointeraction_z <- readRDS("model-fits/z-by-step/growth-models-1-nointeraction-formula.rds")
M1_full_z <- readRDS("model-fits/z-by-step/growth-models-1-full-formula.rds")

M2_nogroup_r <- readRDS("model-fits/raw/growth-models-2-nogroup-formula.rds")
M2_nointeraction_r <- readRDS("model-fits/raw/growth-models-2-nointeraction-formula.rds")
M2_full_r <- readRDS("model-fits/raw/growth-models-2-full-formula.rds")


M2_nogroup_Z <- readRDS("model-fits/z-overall/growth-models-2-nogroup-formula.rds")
M2_nointeraction_Z <- readRDS("model-fits/z-overall/growth-models-2-nointeraction-formula.rds")
M2_full_Z <- readRDS("model-fits/z-overall/growth-models-2-full-formula.rds")

M2_nogroup_z <- readRDS("model-fits/z-by-step/growth-models-2-nogroup-formula.rds")
M2_nointeraction_z <- readRDS("model-fits/z-by-step/growth-models-2-nointeraction-formula.rds")
M2_full_z <- readRDS("model-fits/z-by-step/growth-models-2-full-formula.rds")


# Growth value over Growth Value contrasts ----

## PAC over LOA ----
contrast_pacloa_assoc <- bind_rows(
    nogroup = model_comparison(M2_nogroup_Z$assoc$loa_pac, M1_nogroup_Z$assoc$loa),
    nointeraction = model_comparison(M2_nointeraction_Z$assoc$loa_pac, M1_nointeraction_Z$assoc$loa),
    full = model_comparison(M2_full_Z$assoc$loa_pac, M1_full_Z$assoc$loa),
    .id = "group_cfg"
) |> mutate(contrast = "pac-loa", network = "assoc", scale = "z-overall") |> relocate(contrast, network, scale)

contrast_pacloa_childes <- bind_rows(
    nogroup = model_comparison(M2_nogroup_Z$childes$loa_pac, M1_nogroup_Z$childes$loa),
    nointeraction = model_comparison(M2_nointeraction_Z$childes$loa_pac, M1_nointeraction_Z$childes$loa),
    full = model_comparison(M2_full_Z$childes$loa_pac, M1_full_Z$childes$loa),
    .id = "group_cfg"
) |> mutate(contrast = "pac-loa", network = "childes", scale = "z-overall") |> relocate(contrast, network, scale)

contrast_pacloa_both <- bind_rows(
    nogroup = model_comparison(M2_nogroup_Z$both$loa_pac, M1_nogroup_Z$both$loa),
    nointeraction = model_comparison(M2_nointeraction_Z$both$loa_pac, M1_nointeraction_Z$both$loa),
    full = model_comparison(M2_full_Z$both$loa_pac, M1_full_Z$both$loa),
    .id = "group_cfg"
) |> mutate(contrast = "pac-loa", network = "both", scale = "z-overall") |> relocate(contrast, network, scale)

## LOA over PAC ----
contrast_loapac_assoc <- bind_rows(
    nogroup = model_comparison(M2_nogroup_Z$assoc$loa_pac, M1_nogroup_Z$assoc$pac),
    nointeraction = model_comparison(M2_nointeraction_Z$assoc$loa_pac, M1_nointeraction_Z$assoc$pac),
    full = model_comparison(M2_full_Z$assoc$loa_pac, M1_full_Z$assoc$pac),
    .id = "group_cfg"
) |> mutate(contrast = "loa-pac", network = "assoc", scale = "z-overall") |> relocate(contrast, network, scale)

contrast_loapac_childes <- bind_rows(
    nogroup = model_comparison(M2_nogroup_Z$childes$loa_pac, M1_nogroup_Z$childes$pac),
    nointeraction = model_comparison(M2_nointeraction_Z$childes$loa_pac, M1_nointeraction_Z$childes$pac),
    full = model_comparison(M2_full_Z$childes$loa_pac, M1_full_Z$childes$pac),
    .id = "group_cfg"
) |> mutate(contrast = "loa-pac", network = "childes", scale = "z-overall") |> relocate(contrast, network, scale)

contrast_loapac_both <- bind_rows(
    nogroup = model_comparison(M2_nogroup_Z$both$loa_pac, M1_nogroup_Z$both$pac),
    nointeraction = model_comparison(M2_nointeraction_Z$both$loa_pac, M1_nointeraction_Z$both$pac),
    full = model_comparison(M2_full_Z$both$loa_pac, M1_full_Z$both$pac),
    .id = "group_cfg"
) |> mutate(contrast = "loa-pac", network = "both", scale = "z-overall") |> relocate(contrast, network, scale)

bind_rows(
    contrast_pacloa_assoc, contrast_pacloa_childes, contrast_pacloa_both,
    contrast_loapac_assoc, contrast_loapac_childes, contrast_loapac_both
) |> readr::write_csv("gv_gv_interaction.csv")



model_comparison(M1_full_Z$both$pac, M1_baseinteraction_Z$both$pac)
model_comparison(M1_full_Z$both$loa, M1_baseinteraction_Z$both$loa)

model_comparison(M1_full_Z$assoc$pac, M1_baseinteraction_Z$assoc$pac)
model_comparison(M1_full_Z$assoc$loa, M1_baseinteraction_Z$assoc$loa)

test_all_coefs(M1_full_Z$assoc$pac)
