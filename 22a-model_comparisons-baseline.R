library(dplyr)
library(purrr)
library(tidyr)
library(netgrowr)


ME <- readRDS("model-fits/empty-model.rds")
M0_justgroup <- readRDS("model-fits/baseline-model-justgroup.rds")
M0_nogroup <- readRDS("model-fits/baseline-model-nogroup.rds")
M0_nointeraction <- readRDS("model-fits/baseline-model-nointeraction.rds")
M0_full <- readRDS("model-fits/baseline-model.rds")
MRC_1_nogroup <- readRDS("model-fits/RC-1-model-nogroup.rds")
MRC_1_nointeraction <- readRDS("model-fits/RC-1-model-nointeraction.rds")
MRC_1_full <- readRDS("model-fits/RC-1-model-full.rds")


# Does adding group improve model fit? ----
contrast_group <- bind_rows(
    full = model_comparison(M0_justgroup, ME),
    .id = "group_cfg"
) |> mutate(contrast = "group-empty") |> relocate(contrast)

contrast_baseline <- bind_rows(
    nogroup = model_comparison(M0_nogroup, ME),
    .id = "group_cfg"
) |> mutate(contrast = "B-E") |> relocate(contrast)

contrast_baselinegroup <- bind_rows(
    nointeraction = model_comparison(M0_nointeraction, M0_justgroup),
    .id = "group_cfg"
) |> mutate(contrast = "B+G-G") |> relocate(contrast)

contrast_baselineinteraction <- bind_rows(
    full = model_comparison(M0_full, M0_nointeraction),
    .id = "group_cfg"
) |> mutate(contrast = "BxG-B+G") |> relocate(contrast)

contrast_frequency <- bind_rows(
    nogroup = model_comparison(MRC_1_nogroup$RC1, ME),
    nointeraction = model_comparison(MRC_1_nointeraction$RC1, ME),
    full = model_comparison(MRC_1_full$RC1, ME),
    .id = "group_cfg"
) |> mutate(contrast = "frequency-empty") |> relocate(contrast)

contrast_frequencygroup <- bind_rows(
    nointeraction = model_comparison(MRC_1_nointeraction$RC1, M0_justgroup),
    full = model_comparison(MRC_1_full$RC1, M0_justgroup),
    .id = "group_cfg"
) |> mutate(contrast = "frequency-group") |> relocate(contrast)

contrast_biphon <- bind_rows(
    nogroup = model_comparison(MRC_1_nogroup$RC2, ME),
    nointeraction = model_comparison(MRC_1_nointeraction$RC2, ME),
    full = model_comparison(MRC_1_full$RC2, ME),
    .id = "group_cfg"
) |> mutate(contrast = "biphon-empty") |> relocate(contrast)

contrast_biphongroup <- bind_rows(
    nointeraction = model_comparison(MRC_1_nointeraction$RC2, M0_justgroup),
    full = model_comparison(MRC_1_full$RC2, M0_justgroup),
    .id = "group_cfg"
) |> mutate(contrast = "biphon-group") |> relocate(contrast)

contrast_pndc <- bind_rows(
    nogroup = model_comparison(MRC_1_nogroup$RC3, ME),
    nointeraction = model_comparison(MRC_1_nointeraction$RC3, ME),
    full = model_comparison(MRC_1_full$RC3, ME),
    .id = "group_cfg"
) |> mutate(contrast = "pndc-empty") |> relocate(contrast)

contrast_pndcgroup <- bind_rows(
    nointeraction = model_comparison(MRC_1_nointeraction$RC3, M0_justgroup),
    full = model_comparison(MRC_1_full$RC3, M0_justgroup),
    .id = "group_cfg"
) |> mutate(contrast = "pndc-empty") |> relocate(contrast)
