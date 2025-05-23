library(dplyr)
library(purrr)
library(tidyr)
library(netgrowr)

# Helper functions ----
gvals_over_baseline <- function(gvals, baseline) {
    map2(gvals, baseline, \(m_full_list, m_reduced) {
        map_depth(m_full_list, 2, \(m_full, m_redux) {
            model_comparison(m_full, m_redux) |>
                as.list() |>
                as_tibble()
        }, m_redux = m_reduced) |>
            map(~ list_rbind(.x, names_to = "contrast")) |>
            list_rbind(names_to = "network")
    }) |>
        list_rbind(names_to = "group") |>
        mutate(contrast = paste0(contrast, ">base"))
}

gvals_over_gvals <- function(full, reduced) {
    map2(full, reduced, ~ {
        map2(.x, .y, ~ {
            bind_rows(
                `pac>loa` = model_comparison(.x$loa_pac, .y$loa),
                `loa>pac` = model_comparison(.x$loa_pac, .y$pac),
                .id = "contrast"
            )
        }) |> list_rbind(names_to = "network")
    }) |> list_rbind(names_to = "group")
}

ME <- readRDS("model-fits-split/empty-model.rds")
M0 <- readRDS("model-fits-split/baseline-model.rds")
M1_r <- readRDS("model-fits-split/raw/growth-model-1.rds")
M2_r <- readRDS("model-fits-split/raw/growth-model-2.rds")
M3_r <- readRDS("model-fits-split/raw/growth-model-3.rds")
M1_Z <- readRDS("model-fits-split/z-overall/growth-model-1.rds")
M2_Z <- readRDS("model-fits-split/z-overall/growth-model-2.rds")
M3_Z <- readRDS("model-fits-split/z-overall/growth-model-3.rds")
M1_z <- readRDS("model-fits-split/z-by-step/growth-model-1.rds")
M2_z <- readRDS("model-fits-split/z-by-step/growth-model-2.rds")
M3_z <- readRDS("model-fits-split/z-by-step/growth-model-3.rds")

# Baseline over empty ----
df_base_empty <- map2(M0, ME, ~ {
    model_comparison(.x, .y) |>
        as.list() |>
        as_tibble()
}) |> list_rbind(names_to = "group") |>
    mutate(contrast = "base>empty") |>
    relocate(contrast, .after = group)


# 1 GV over baseline ----
df_gv_base <- bind_rows(
    raw = gvals_over_baseline(M1_r, M0),
    z_overall = gvals_over_baseline(M1_Z, M0),
    z_by_step = gvals_over_baseline(M1_z, M0),
    .id = "scale"
)

df_gv_base |>
    filter(
        network == "both",
        scale == "z_overall"
    )

df_gv_base |>
    filter(
        p < .05,
        network == "both",
        scale == "z_by_step"
    )


# GV over GV ----
df_gv_gv <- bind_rows(
    raw = gvals_over_gvals(M2_r, M1_r),
    z_overall = gvals_over_gvals(M2_Z, M1_Z),
    z_by_step = gvals_over_gvals(M2_z, M1_z),
    .id = "scale"
)


df_gv_gv |>
    filter(
        p < .05,
        network == "both",
        scale == "z_overall"
    )

df_gv_gv |>
    filter(
        p < .05,
        network == "both",
        scale == "z_by_step"
    )

df_gv_gv |>
    filter(
        network == "both",
        scale == "raw"
    )
