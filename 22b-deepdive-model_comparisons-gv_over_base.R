library(dplyr)
library(purrr)
library(tidyr)
library(netgrowr)


M0_nogroup <- readRDS("model-fits/baseline-model-nogroup.rds")
M0_nointeraction <- readRDS("model-fits/baseline-model-nointeraction.rds")
M0_full <- readRDS("model-fits/baseline-model.rds")

M1_nogroup_r <- readRDS("model-fits/raw/growth-models-1-nogroup.rds")
M1_nointeraction_r <- readRDS("model-fits/raw/growth-models-1-nointeraction.rds")
M1_baseinteraction_r <- readRDS("model-fits/raw/growth-models-1-baseinteraction.rds")
M1_full_r <- readRDS("model-fits/raw/growth-models-1-full.rds")

M1_nogroup_Z <- readRDS("model-fits/z-overall/growth-models-1-nogroup.rds")
M1_nointeraction_Z <- readRDS("model-fits/z-overall/growth-models-1-nointeraction.rds")
M1_baseinteraction_Z <- readRDS("model-fits/z-overall/growth-models-1-baseinteraction.rds")
M1_full_Z <- readRDS("model-fits/z-overall/growth-models-1-full.rds")

M1_nogroup_z <- readRDS("model-fits/z-by-step/growth-models-1-nogroup.rds")
M1_nointeraction_z <- readRDS("model-fits/z-by-step/growth-models-1-nointeraction.rds")
M1_baseinteraction_z <- readRDS("model-fits/z-by-step/growth-models-1-baseinteraction.rds")
M1_full_z <- readRDS("model-fits/z-by-step/growth-models-1-full.rds")


# Growth Value over Baseline contrasts ----
contrast_patbaseline_assoc <- bind_rows(
    nogroup = model_comparison(M1_nogroup_z$assoc$pat, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_z$assoc$pat, M0_nointeraction),
    full = model_comparison(M1_full_z$assoc$pat, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "pat-baseline", network = "assoc", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_patbaseline_childes <- bind_rows(
    nogroup = model_comparison(M1_nogroup_z$childes$pat, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_z$childes$pat, M0_nointeraction),
    full = model_comparison(M1_full_z$childes$pat, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "pat-baseline", network = "childes", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_patbaseline_both <- bind_rows(
    nogroup = model_comparison(M1_nogroup_z$both$pat, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_z$both$pat, M0_nointeraction),
    full = model_comparison(M1_full_z$both$pat, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "pat-baseline", network = "both", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_loabaseline_assoc <- bind_rows(
    nogroup = model_comparison(M1_nogroup_z$assoc$loa, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_z$assoc$loa, M0_nointeraction),
    full = model_comparison(M1_full_z$assoc$loa, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "loa-baseline", network = "assoc", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_loabaseline_childes <- bind_rows(
    nogroup = model_comparison(M1_nogroup_z$childes$loa, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_z$childes$loa, M0_nointeraction),
    full = model_comparison(M1_full_z$childes$loa, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "loa-baseline", network = "childes", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_loabaseline_both <- bind_rows(
    nogroup = model_comparison(M1_nogroup_z$both$loa, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_z$both$loa, M0_nointeraction),
    full = model_comparison(M1_full_z$both$loa, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "loa-baseline", network = "both", scale = "z-by-step") |> relocate(contrast, network, scale)


contrast_pacbaseline_assoc <- bind_rows(
    nogroup = model_comparison(M1_nogroup_z$assoc$pac, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_z$assoc$pac, M0_nointeraction),
    full = model_comparison(M1_full_z$assoc$pac, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "pac-baseline", network = "assoc", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_pacbaseline_childes <- bind_rows(
    nogroup = model_comparison(M1_nogroup_z$childes$pac, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_z$childes$pac, M0_nointeraction),
    full = model_comparison(M1_full_z$childes$pac, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "pac-baseline", network = "childes", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_pacbaseline_both <- bind_rows(
    nogroup = model_comparison(M1_nogroup_z$both$pac, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_z$both$pac, M0_nointeraction),
    full = model_comparison(M1_full_z$both$pac, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "pac-baseline", network = "both", scale = "z-by-step") |> relocate(contrast, network, scale)

bind_rows(
    contrast_patbaseline_assoc, contrast_patbaseline_childes, contrast_patbaseline_both,
    contrast_loabaseline_assoc, contrast_loabaseline_childes, contrast_loabaseline_both,
    contrast_pacbaseline_assoc, contrast_pacbaseline_childes, contrast_pacbaseline_both,
)

model_comparison(M1_baseinteraction_z$assoc$pac, M0_full) |> as.list() |> as_tibble()









# Growth Value over Baseline contrasts (z-by-step) ----
contrast_patbaseline_assoc <- bind_rows(
    nogroup = model_comparison(M1_nogroup_Z$assoc$pat, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_Z$assoc$pat, M0_nointeraction),
    full = model_comparison(M1_full_Z$assoc$pat, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "pat-baseline", network = "assoc", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_patbaseline_childes <- bind_rows(
    nogroup = model_comparison(M1_nogroup_Z$childes$pat, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_Z$childes$pat, M0_nointeraction),
    full = model_comparison(M1_full_Z$childes$pat, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "pat-baseline", network = "childes", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_patbaseline_both <- bind_rows(
    nogroup = model_comparison(M1_nogroup_Z$both$pat, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_Z$both$pat, M0_nointeraction),
    full = model_comparison(M1_full_Z$both$pat, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "pat-baseline", network = "both", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_loabaseline_assoc <- bind_rows(
    nogroup = model_comparison(M1_nogroup_Z$assoc$loa, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_Z$assoc$loa, M0_nointeraction),
    full = model_comparison(M1_full_Z$assoc$loa, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "loa-baseline", network = "assoc", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_loabaseline_childes <- bind_rows(
    nogroup = model_comparison(M1_nogroup_Z$childes$loa, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_Z$childes$loa, M0_nointeraction),
    full = model_comparison(M1_full_Z$childes$loa, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "loa-baseline", network = "childes", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_loabaseline_both <- bind_rows(
    nogroup = model_comparison(M1_nogroup_Z$both$loa, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_Z$both$loa, M0_nointeraction),
    full = model_comparison(M1_full_Z$both$loa, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "loa-baseline", network = "both", scale = "z-by-step") |> relocate(contrast, network, scale)


contrast_pacbaseline_assoc <- bind_rows(
    nogroup = model_comparison(M1_nogroup_Z$assoc$pac, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_Z$assoc$pac, M0_nointeraction),
    full = model_comparison(M1_full_Z$assoc$pac, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "pac-baseline", network = "assoc", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_pacbaseline_childes <- bind_rows(
    nogroup = model_comparison(M1_nogroup_Z$childes$pac, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_Z$childes$pac, M0_nointeraction),
    full = model_comparison(M1_full_Z$childes$pac, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "pac-baseline", network = "childes", scale = "z-by-step") |> relocate(contrast, network, scale)

contrast_pacbaseline_both <- bind_rows(
    nogroup = model_comparison(M1_nogroup_Z$both$pac, M0_nogroup),
    nointeraction = model_comparison(M1_nointeraction_Z$both$pac, M0_nointeraction),
    full = model_comparison(M1_full_Z$both$pac, M0_full),
    .id = "group_cfg"
) |> mutate(contrast = "pac-baseline", network = "both", scale = "z-by-step") |> relocate(contrast, network, scale)

bind_rows(
    contrast_patbaseline_assoc, contrast_patbaseline_childes, contrast_patbaseline_both,
    contrast_loabaseline_assoc, contrast_loabaseline_childes, contrast_loabaseline_both,
    contrast_pacbaseline_assoc, contrast_pacbaseline_childes, contrast_pacbaseline_both,
) |> readr::write_csv("gv_baseline_interaction.csv")

model_comparison(M1_baseinteraction_Z$assoc$pac, M0_full) |> as.list() |> as_tibble()







