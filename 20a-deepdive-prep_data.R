library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)


# Rename variables ----
modelvars_df <- readRDS("./network/modelvars_vsoa_2025_05_14_relearn0.rds") |>
    mutate(
        label = paste(group, num_item_id)
    ) |>
    rename(
        vsoa = aoa,
        step = month,
        gv_pat1 = preferential_attachment_child,
        gv_pat2 = preferential_attachment_childes,
        gv_pac1 = preferential_acquisition_child,
        gv_pac2 = preferential_acquisition_childes,
        gv_loa1 = lure_of_the_associates_child,
        gv_loa2 = lure_of_the_associates_childes
    )


# Standardize growth values by step ----
modelvars_df <- modelvars_df |>
    group_by(group, step) |>
    mutate(
        across(starts_with("gv_"), list(z = ~ (.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)))
    )

modelvars_df |>
    select(group, num_item_id, word, step, vsoa, learned, starts_with("gv_")) |>
    pivot_longer(
        cols = starts_with("gv_"),
        names_to = c("growth_model", "source"),
        names_pattern = "gv_([pl][ao][act])([12])",
        values_to = "growth_value"
    ) |>
    mutate(source = factor(as.numeric(source), 1:2, c("assoc", "childes"))) |>
    drop_na() |>
    group_by(group, step, learned, growth_model, source) |>
    summarize(
        across(growth_value, list(nnz = ~{sum(.x != 0)}, sum = sum, mean = mean, sd = sd, sumz = ~{sum(.x) / sd(.x)}, max = max))
    ) |>
    filter(learned == FALSE) |>
    ggplot(aes(x = step, y = growth_value_max, color = growth_model)) +
        geom_point() +
        facet_grid(group ~ source)


confounds <- modelvars_df |>
    select(num_item_id, word, nphon, CHILDES_Freq, BiphonProb.avg, PNDC.avg) |>
    distinct() |>
    mutate(across(c(nphon, CHILDES_Freq), list(log = ~ log(.x + 1))))


# Orthogonalizing Confounds ----

## correlation ----
# Number of phonemes is strongly (anti)correlated with frequency and
# phonological neighborhood density. It is moderately correlated with average
# biphone probability. Additionally, neighborhood density is correlated with frequency.
confounds |>
    select(nphon_log, CHILDES_Freq_log, BiphonProb.avg, PNDC.avg) |>
    cor()

## pca ----
# A principal components analysis which selects and rotates the first three
# factors obtains a solution where, essentially, number of phonemes is folded
# into each factor proportional to its correlation, and the three factors
# correspond more or less with Frequency, neighborhood density, and biphon
# probability. These rotated components may be a little cleaner to work with.
pca_varimax <- principal(
    select(confounds, nphon_log, CHILDES_Freq_log, BiphonProb.avg, PNDC.avg),
    nfactors = 3,
    rotate = "varimax"
)


## update dataframe ----
confounds <- bind_cols(
    confounds,
    pca_varimax$scores |> as_tibble()
)
modelvars_df <- modelvars_df |>
    left_join(select(confounds, num_item_id, RC1, RC2, RC3)) |>
    select(-nphon, -CHILDES_Freq, -BiphonProb.avg, -PNDC.avg) |>
    relocate(label, num_item_id, .before = word) |>
    rename_with(
        ~ stringr::str_remove(.x, "gv_"),
        starts_with("gv_")
    )


# Set contrast code for group variable ----
contrasts(modelvars_df$group) <- c(-.5, .5)


# Save revised dataframe ----
saveRDS(modelvars_df, "network/modelvars_vsoa_RC_z.rds")


# Split dataframe by group ----
modelvars <- modelvars_df |>
    group_by(group) |>
    group_split()

names(modelvars) <- modelvars_df |>
    group_by(group) |>
    group_keys() |>
    pull(group)


## Append group label to growth model variables ----
modelvars_rn <- imap(modelvars, function(.x, group) {
    rename_with(
        .data = .x,
        .fn = function(x, y) {
            paste(x, y, sep = "_")
        },
        .cols = c(
            pat1,
            pat2,
            loa1,
            loa2,
            pac1,
            pac2
        ),
        y = group
    )
})


## standardize growth values ----
modelvars_rn <- imap(modelvars_rn, function(x, g) {
    x |>
        mutate(across(ends_with(g), list(z_overall = ~ (.x - mean(.x, na.rm = T)) / sd(.x, na.rm = T)))) |>
        group_by(month) |>
        mutate(across(ends_with(g), list(z = ~ (.x - mean(.x, na.rm = T)) / sd(.x, na.rm = T))))
})

saveRDS(modelvars_rn, "network/modelvars_vsoa_RC_z_split.rds")
