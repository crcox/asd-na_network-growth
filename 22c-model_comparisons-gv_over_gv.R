library(dplyr)
library(purrr)
library(tidyr)
library(netgrowr)


M1_nogroup_r <- readRDS("model-fits/raw/growth-models-1-nogroup.rds")
M1_nointeraction_r <- readRDS("model-fits/raw/growth-models-1-nointeraction.rds")
M1_full_r <- readRDS("model-fits/raw/growth-models-1-full.rds")

M1_nogroup_Z <- readRDS("model-fits/z-overall/growth-models-1-nogroup.rds")
M1_nointeraction_Z <- readRDS("model-fits/z-overall/growth-models-1-nointeraction.rds")
M1_full_Z <- readRDS("model-fits/z-overall/growth-models-1-full.rds")

M1_nogroup_z <- readRDS("model-fits/z-by-step/growth-models-1-nogroup.rds")
M1_nointeraction_z <- readRDS("model-fits/z-by-step/growth-models-1-nointeraction.rds")
M1_full_z <- readRDS("model-fits/z-by-step/growth-models-1-full.rds")

M2_nogroup_r <- readRDS("model-fits/raw/growth-models-2-nogroup.rds")
M2_nointeraction_r <- readRDS("model-fits/raw/growth-models-2-nointeraction.rds")
M2_full_r <- readRDS("model-fits/raw/growth-models-2-full.rds")

M2_nogroup_Z <- readRDS("model-fits/z-overall/growth-models-2-nogroup.rds")
M2_nointeraction_Z <- readRDS("model-fits/z-overall/growth-models-2-nointeraction.rds")
M2_full_Z <- readRDS("model-fits/z-overall/growth-models-2-full.rds")

M2_nogroup_z <- readRDS("model-fits/z-by-step/growth-models-2-nogroup.rds")
M2_nointeraction_z <- readRDS("model-fits/z-by-step/growth-models-2-nointeraction.rds")
M2_full_z <- readRDS("model-fits/z-by-step/growth-models-2-full.rds")


# Growth value over Growth Value contrasts ----

## PAC over LOA ----
contrast_pacloa_assoc <- bind_rows(
    nogroup = model_comparison(M2_nogroup_z$assoc$loa_pac, M1_nogroup_z$assoc$loa),
    nointeraction = model_comparison(M2_nointeraction_z$assoc$loa_pac, M1_nointeraction_z$assoc$loa),
    full = model_comparison(M2_full_z$assoc$loa_pac, M1_full_z$assoc$loa),
    .id = "group_cfg"
) |> mutate(contrast = "pac-loa", network = "assoc", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_pacloa_childes <- bind_rows(
    nogroup = model_comparison(M2_nogroup_z$childes$loa_pac, M1_nogroup_z$childes$loa),
    nointeraction = model_comparison(M2_nointeraction_z$childes$loa_pac, M1_nointeraction_z$childes$loa),
    full = model_comparison(M2_full_z$childes$loa_pac, M1_full_z$childes$loa),
    .id = "group_cfg"
) |> mutate(contrast = "pac-loa", network = "childes", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_pacloa_both <- bind_rows(
    nogroup = model_comparison(M2_nogroup_z$both$loa_pac, M1_nogroup_z$both$loa),
    nointeraction = model_comparison(M2_nointeraction_z$both$loa_pac, M1_nointeraction_z$both$loa),
    full = model_comparison(M2_full_z$both$loa_pac, M1_full_z$both$loa),
    .id = "group_cfg"
) |> mutate(contrast = "pac-loa", network = "both", scale = "z-by-step") |> relocate(contrast, network, scale)

## LOA over PAC ----
contrast_loapac_assoc <- bind_rows(
    nogroup = model_comparison(M2_nogroup_z$assoc$loa_pac, M1_nogroup_z$assoc$pac),
    nointeraction = model_comparison(M2_nointeraction_z$assoc$loa_pac, M1_nointeraction_z$assoc$pac),
    full = model_comparison(M2_full_z$assoc$loa_pac, M1_full_z$assoc$pac),
    .id = "group_cfg"
) |> mutate(contrast = "loa-pac", network = "assoc", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_loapac_childes <- bind_rows(
    nogroup = model_comparison(M2_nogroup_z$childes$loa_pac, M1_nogroup_z$childes$pac),
    nointeraction = model_comparison(M2_nointeraction_z$childes$loa_pac, M1_nointeraction_z$childes$pac),
    full = model_comparison(M2_full_z$childes$loa_pac, M1_full_z$childes$pac),
    .id = "group_cfg"
) |> mutate(contrast = "loa-pac", network = "childes", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_loapac_both <- bind_rows(
    nogroup = model_comparison(M2_nogroup_z$both$loa_pac, M1_nogroup_z$both$pac),
    nointeraction = model_comparison(M2_nointeraction_z$both$loa_pac, M1_nointeraction_z$both$pac),
    full = model_comparison(M2_full_z$both$loa_pac, M1_full_z$both$pac),
    .id = "group_cfg"
) |> mutate(contrast = "loa-pac", network = "both", scale = "z-by-step") |> relocate(contrast, network, scale)

bind_rows(
    contrast_pacloa_assoc, contrast_pacloa_childes, contrast_pacloa_both,
    contrast_loapac_assoc, contrast_loapac_childes, contrast_loapac_both
)
