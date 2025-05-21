library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(netgrowr)
library(parallel)
library(progressr)

source("./R/utils.R")

append_formula <- function(f, x) {
    update(f, paste("~.", x, sep = "+"))
}

combn_and_paste <- function(x, m) {
    return(apply(combn(x, m), MARGIN = 2, paste, collapse = " + "))
}

model_comp_helper <- function(full, restricted, M) {
    return(netgrowr::model_comparison(M[[full]], M[[restricted]]))
}

load_variable <- function(file, variable) {
    e <- new.env()
    load(file, envir = e)
    return(e[[variable]])
}

list_or_vec <- function(x) {
    return(is.list(x) || is.vector(x))
}

# Load data ----
# Month 30 is the highest month we have AoA data for.
# If we model month 30, all unknown words will be learned.
# To avoid this, month 30 is not modeled.
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
    ) |>
    group_by(group, step) |>
    mutate(
        across(starts_with("gv_"), list(z = ~ (.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)))
    )

modelvars_df |>
    select(group, vs, num_item_id, word, month, aoa, learned, ends_with("_child"), ends_with("_childes")) |>
    rename_with(~{sub("_ch", "Qch", .x)}) |>
    pivot_longer(
        cols = contains("Q"),
        names_to = c("growth_model", "source"),
        names_sep = "Q",
        values_to = "growth_value"
    ) |>
    drop_na(growth_value) |>
    group_by(group, vs, month, learned, growth_model, source) |>
    summarize(
        across(growth_value, list(nnz = ~{sum(.x != 0)}, sum = sum, mean = mean, sd = sd, sumz = ~{sum(.x) / sd(.x)}, max = max))
    ) |>
    group_by(group, learned, growth_model, source) |>
    mutate(growth_value_sum_scaled = growth_value_sum / max(growth_value_sum)) |>
    ungroup() |>
    filter(learned == FALSE) |>
    ggplot(aes(x = vs, y = growth_value_sum_scaled, color = growth_model)) +
        geom_point() +
        facet_grid(group ~ source) +
        xlab("vocabulary size") +
        ylab("scaled sum of growth values\nover unknown words") +
        theme_bw(base_size = 14) +
        scale_color_discrete(labels = c("LOA", "PAC", "PAT")) +
        labs(color = "growth model")


modelvars <- modelvars_df |>
    group_by(group) |>
    group_split()

names(modelvars) <- modelvars_df |>
    group_by(group) |>
    group_keys() |>
    pull(group)


# Append group label to growth model variables
modelvars_rn <- imap(modelvars, function(.x, group) {
    rename_with(
        .data = .x,
        .fn = function(x, y) {
            paste(x, y, sep = "_")
        },
        .cols = c(
            preferential_attachment_child,
            lure_of_the_associates_child,
            preferential_acquisition_child,
            preferential_attachment_childes,
            lure_of_the_associates_childes,
            preferential_acquisition_childes
        ),
        y = group
    )
})

modelvars_rn <- imap(modelvars_rn, function(x, g) {
    x |>
        mutate(across(ends_with(g), list(z = Z))) |>
        group_by(month) |>
        mutate(across(ends_with(g), list(z_month = Z)))
})

confounds <- modelvars_df |>
    select(num_item_id, word, nphon, CHILDES_Freq, BiphonProb.avg, PNDC.avg) |>
    distinct() |>
    mutate(across(c(nphon, CHILDES_Freq), list(log = ~ log(.x + 1))))

# Orthogonalizing Confounds ----
# Number of phonemes is strongly (anti)correlated with frequency and
# phonological neighborhood density. It is moderately correlated with average
# biphone probability. Additionally, neighborhood density is correlated with frequency.
confounds |>
    select(nphon_log, CHILDES_Freq_log, BiphonProb.avg, PNDC.avg) |>
    cor()

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



# Define models by constructing formulas ----

# Empty (random selection) model
# +++ All words assigned equal probability
fE <- formula(vsoa ~ 1)
ME <- netgrowr::mle_network_growth(fE, data = na.omit(modelvars_df), split_by = "month", label_with = "label")

# Psycholinguistic baseline model
# +++ number of phonemes
# +++ CHILDES Frequency
# +++ phonotactic probability (biphone)
# +++ phonological neighborhood density (KU child corpus)
f0 <- update(fE, ~ . + RC1 + RC2 + RC3)
M0_split <- map(modelvars, ~ {
    netgrowr::mle_network_growth(f0, data = na.omit(.x), split_by = "month", label_with = "num_item_id")
})

M0 <- netgrowr::mle_network_growth(
    f0,
    data = modelvars_df,
    split_by = "month",
    label_with = "label"
)

contrasts(modelvars_df$group) <- c(-.5, .5)

f0group <- update(f0, ~ . + group)
M0group <- netgrowr::mle_network_growth(
    f0group,
    data = modelvars_df,
    split_by = "month",
    label_with = "label"
)

f0groupx <- update(f0, ~ (.) * group)
M0groupx <- netgrowr::mle_network_growth(
    f0groupx,
    data = modelvars_df,
    split_by = "month",
    label_with = "label"
)

f1group <- list(
    assoc = list(
        pat = update(f0, ~ . + pat1_z + group),
        loa = update(f0, ~ . + loa1_z + group),
        pac = update(f0, ~ . + pac1_z + group)
    ),
    childes = list(
        pat = update(f0, ~ . + pat2_z + group),
        loa = update(f0, ~ . + loa2_z + group),
        pac = update(f0, ~ . + pac2_z + group)
    ),
    both = list(
        pat = update(f0, ~ . + pat1_z + pat2_z + group ),
        loa = update(f0, ~ . + loa1_z + loa2_z + group ),
        pac = update(f0, ~ . + pac1_z + pac2_z + group )
    )
)

M1group <- map_depth(f1group, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "step",
        label_with = "label"
    )

}, .progress = TRUE)

model_comparison(M1group$assoc, M0group)


f1groupx <- list(
    assoc = list(
        pat = update(f0, ~ (. + pat1_z) * group),
        loa = update(f0, ~ (. + loa1_z) * group),
        pac = update(f0, ~ (. + pac1_z) * group)
    ),
    childes = list(
        pat = update(f0, ~ (. + pat2_z) * group),
        loa = update(f0, ~ (. + loa2_z) * group),
        pac = update(f0, ~ (. + pac2_z) * group)
    ),
    both = list(
        pat = update(f0, ~ (. + pat1_z + pat2_z) * group ),
        loa = update(f0, ~ (. + loa1_z + loa2_z) * group ),
        pac = update(f0, ~ (. + pac1_z + pac2_z) * group )
    )
)
M1groupx <- map_depth(f1groupx, 2, ~{
    netgrowr::mle_network_growth(
        .x,
        data = modelvars_df,
        split_by = "step",
        label_with = "label"
    )

}, .progress = TRUE)


f0_nofreq <- update(fE, ~ . + I(Z(log(nphon + 1))) + Z(BiphonProb.avg) + Z(PNDC.avg))
M0_nofreq <- map(modelvars, ~ {
    netgrowr::mle_network_growth(f0_nofreq, data = na.omit(.x), split_by = "month", label_with = "num_item_id")
})

# Networks over control
network_gv1 <- list(
    autistic = list(
        assoc = c(
            pat = "preferential_attachment_child_autistic",
            pac = "preferential_acquisition_child_autistic",
            loa = "lure_of_the_associates_child_autistic"
        ),
        childes = c(
            pat = "preferential_attachment_childes_autistic",
            pac = "preferential_acquisition_childes_autistic",
            loa = "lure_of_the_associates_childes_autistic"
        ),
        both = c(
            pat = "preferential_attachment_child_autistic + preferential_attachment_childes_autistic",
            pac = "preferential_acquisition_child_autistic + preferential_acquisition_childes_autistic",
            loa = "lure_of_the_associates_child_autistic + lure_of_the_associates_childes_autistic"
        )
    ),
    nonautistic = list(
        assoc = c(
            pat = "preferential_attachment_child_nonautistic",
            pac = "preferential_acquisition_child_nonautistic",
            loa = "lure_of_the_associates_child_nonautistic"
        ),
        childes = c(
            pat = "preferential_attachment_childes_nonautistic",
            pac = "preferential_acquisition_childes_nonautistic",
            loa = "lure_of_the_associates_childes_nonautistic"
        ),
        both = c(
            pat = "preferential_attachment_child_nonautistic + preferential_attachment_childes_nonautistic",
            pac = "preferential_acquisition_child_nonautistic + preferential_acquisition_childes_nonautistic",
            loa = "lure_of_the_associates_child_nonautistic + lure_of_the_associates_childes_nonautistic"
        )
    )
)
network_gv1 <- list(
    autistic = list(
        assoc = c(
            pac = "preferential_acquisition_child_autistic"
        ),
        childes = c(
            pac = "preferential_acquisition_childes_autistic"
        ),
        both = c(
            pac = "preferential_acquisition_child_autistic + preferential_acquisition_childes_autistic"
        )
    ),
    nonautistic = list(
        assoc = c(
            pac = "preferential_acquisition_child_nonautistic"
        ),
        childes = c(
            pac = "preferential_acquisition_childes_nonautistic"
        ),
        both = c(
            pac = "preferential_acquisition_child_nonautistic + preferential_acquisition_childes_nonautistic"
        )
    )
)
network_gv1 <- list(
    autistic = list(
        assoc = c(
            pac = "preferential_acquisition_child_autistic_z"
        ),
        childes = c(
            pac = "preferential_acquisition_childes_autistic_z"
        ),
        both = c(
            pac = "preferential_acquisition_child_autistic_z + preferential_acquisition_childes_autistic_z"
        )
    ),
    nonautistic = list(
        assoc = c(
            pac = "preferential_acquisition_child_nonautistic_z"
        ),
        childes = c(
            pac = "preferential_acquisition_childes_nonautistic_z"
        ),
        both = c(
            pac = "preferential_acquisition_child_nonautistic_z + preferential_acquisition_childes_nonautistic_z"
        )
    )
)
f1 <- map_depth(network_gv1, 3, append_formula, f = f0, .is_node = list_or_vec, .progress = TRUE)


m <- modelvars_rn$nonautistic |>
    group_by(month)

modelvars_resid <- map(group_split(m), \(x) {
    z <- !is.na(x$preferential_acquisition_child_nonautistic)
    x$nphon[z] <- resid(lm(I(Z(log(nphon + 1))) ~ preferential_acquisition_child_nonautistic, data = x))
    x$CHILDES_Freq[z] <- resid(lm(I(Z(log(CHILDES_Freq + 1))) ~ preferential_acquisition_child_nonautistic, data = x))
    x$BiphonProb.avg[z] <- resid(lm(Z(BiphonProb.avg) ~ preferential_acquisition_child_nonautistic, data = x))
    x$PNDC.avg[z] <- resid(lm(Z(PNDC.avg) ~ preferential_acquisition_child_nonautistic, data = x))
    return(x)
}) |>
    list_rbind()

modelvars_resid_childes <- map(group_split(m), \(x) {
    z <- !is.na(x$preferential_acquisition_child_nonautistic)
    x$nphon[z] <- resid(lm(I(Z(log(nphon + 1))) ~ preferential_acquisition_childes_nonautistic, data = x))
    x$CHILDES_Freq[z] <- resid(lm(I(Z(log(CHILDES_Freq + 1))) ~ preferential_acquisition_childes_nonautistic, data = x))
    x$BiphonProb.avg[z] <- resid(lm(Z(BiphonProb.avg) ~ preferential_acquisition_childes_nonautistic, data = x))
    x$PNDC.avg[z] <- resid(lm(Z(PNDC.avg) ~ preferential_acquisition_childes_nonautistic, data = x))
    return(x)
}) |>
    list_rbind()

M0_resid <-  netgrowr::mle_network_growth(aoa ~ nphon + CHILDES_Freq + BiphonProb.avg + PNDC.avg, data = drop_na(modelvars_resid), split_by = "month", label_with = "num_item_id")
M0_resid_childes <-  netgrowr::mle_network_growth(aoa ~ nphon + CHILDES_Freq + BiphonProb.avg + PNDC.avg, data = drop_na(modelvars_resid_childes), split_by = "month", label_with = "num_item_id")

model_comparison(M0_resid, M0$nonautistic)
model_comparison(M0_resid_childes, M0$nonautistic)


M1_z <- map2(modelvars_rn, f1, \(mvars_group, f1_group, p) {
    map_depth(f1_group, .depth = 2, \(f, d, p) {
        p()
        netgrowr::mle_network_growth(f, data = na.omit(d), split_by = "month", label_with = "num_item_id")
    }, d = mvars_group, p = p)
}, p = progressor(6)) |> with_progress()

model_comparison(M1_z$nonautistic$assoc$pac, M0$nonautistic)
model_comparison(M1_z$nonautistic$childes$pac, M0$nonautistic)
model_comparison(M1_z$nonautistic$both$pac, M0$nonautistic)
model_comparison(M1_z$autistic$assoc$pac, M0$autistic)
model_comparison(M1_z$autistic$childes$pac, M0$autistic)
model_comparison(M1_z$autistic$both$pac, M0$autistic)

model_comparison(M1_z$nonautistic$both$pac, M1_z$nonautistic$assoc$pac)
model_comparison(M1_z$nonautistic$both$pac, M1_z$nonautistic$childes$pac)
model_comparison(M1_z$autistic$both$pac, M1_z$autistic$assoc$pac)
model_comparison(M1_z$autistic$both$pac, M1_z$autistic$childes$pac)

M1_z_nofreq <- map2(modelvars_rn, f1, \(mvars_group, f1_group, p) {
    map_depth(f1_group, .depth = 2, \(f, d, p) {
        p()
        netgrowr::mle_network_growth(update(f, "~.-I(Z(log(CHILDES_Freq + 1)))"), data = na.omit(d), split_by = "month", label_with = "num_item_id")
    }, d = mvars_group, p = p)
}, p = progressor(6)) |> with_progress()

model_comparison(M1_z_nofreq$nonautistic$assoc$pac, M0_nofreq$nonautistic)
model_comparison(M1_z_nofreq$nonautistic$childes$pac, M0_nofreq$nonautistic)
model_comparison(M1_z_nofreq$nonautistic$both$pac, M0_nofreq$nonautistic)
model_comparison(M1_z_nofreq$autistic$assoc$pac, M0_nofreq$autistic)
model_comparison(M1_z_nofreq$autistic$childes$pac, M0_nofreq$autistic)
model_comparison(M1_z_nofreq$autistic$both$pac, M0_nofreq$autistic)


model_comparison(M1_z_nofreq$nonautistic$both$pac, M1_z_nofreq$nonautistic$assoc$pac)
model_comparison(M1_z_nofreq$nonautistic$both$pac, M1_z_nofreq$nonautistic$childes$pac)
model_comparison(M1_z_nofreq$autistic$both$pac, M1_z_nofreq$autistic$assoc$pac)
model_comparison(M1_z_nofreq$autistic$both$pac, M1_z_nofreq$autistic$childes$pac)

lab <- map(M1_z, ~ names(.x$both$pac$coefficients))

rbind(
    M0_nofreq$nonautistic$coefficients[lab$nonautistic],
    M1_z_nofreq$nonautistic$assoc$pac$coefficients[lab$nonautistic],
    M1_z_nofreq$nonautistic$childes$pac$coefficients[lab$nonautistic],
    M1_z_nofreq$nonautistic$both$pac$coefficients[lab$nonautistic],
    M0_nofreq$autistic$coefficients[lab$autistic],
    M1_z_nofreq$autistic$assoc$pac$coefficients[lab$autistic],
    M1_z_nofreq$autistic$childes$pac$coefficients[lab$autistic],
    M1_z_nofreq$autistic$both$pac$coefficients[lab$autistic]
) |>unname()

model_comparison(M1$nonautistic$assoc$pac, M0_resid)
model_comparison(M1$nonautistic$childes$pac, M0_resid_childes)
model_comparison(M1$nonautistic$both$pac, M0_resid_childes)

M1 |>
    map_dfr(\(g) {
        map_dfr(g, \(s) {
            map_dfr(s, ~ tibble(regressor = names(coef(.x)), b = coef(.x)), .id = "model")
        }, .id = "source")
    }, .id = "group") |>
    mutate(
        across(c(group, source, model), as.factor),
        regressor = factor(regressor, levels = c(
            "(intercept)",
            "I(Z(log(nphon + 1)))",
            "I(Z(log(CHILDES_Freq + 1)))",
            "Z(BiphonProb.avg)",
            "Z(PNDC.avg)",
            "preferential_attachment_child_autistic",
            "preferential_attachment_childes_autistic",
            "preferential_attachment_child_nonautistic",
            "preferential_attachment_childes_nonautistic",
            "preferential_acquisition_child_autistic",
            "preferential_acquisition_childes_autistic",
            "preferential_acquisition_child_nonautistic",
            "preferential_acquisition_childes_nonautistic",
            "lure_of_the_associates_child_autistic",
            "lure_of_the_associates_childes_autistic",
            "lure_of_the_associates_child_nonautistic",
            "lure_of_the_associates_childes_nonautistic"
        ), labels = c("b0", "nphon", "freq", "biphon prob", "pndc", "pat", "pat", "pat", "pat", "pac", "pac", "pac", "pac", "loa", "loa", "loa", "loa"))
    )

saveRDS(M0, "fits-with-no-growth-models-20250514.rds")
saveRDS(M1, "fits-with-one-growth-model-20250514.rds")
saveRDS(modelvars_rn, "model-vars-20250514.rds")


# Models that combine growth values derived from different growth models for the same network
network_gv2 <- list(
    autistic = combn_and_paste(c(
        "preferential_attachment_child_autistic + preferential_attachment_childes_autistic",
        "preferential_acquisition_child_autistic + preferential_acquisition_childes_autistic",
        "lure_of_the_associates_child_autistic + lure_of_the_associates_childes_autistic"
    ), m = 2),
    nonautistic = combn_and_paste(c(
        "preferential_attachment_child_nonautistic + preferential_attachment_childes_nonautistic",
        "preferential_acquisition_child_nonautistic + preferential_acquisition_childes_nonautistic",
        "lure_of_the_associates_child_nonautistic + lure_of_the_associates_childes_nonautistic"
    ), m = 2)
)
f2g <- map_depth(network_gv2, 2, append_formula, f = f0, .is_node = list_or_vec)

# Define cluster for parallel computing and optimize models ----

# The optimization process is not guaranteed to find a global minimum.
# Different runs of the optimization will yield numeric differences, and
# rarely very bad solutions are obtained. Thus, the optimization process is
# repeated 500 times for every model fit. The minimum negative log likelihood
# for each model over the 500 optimization attempts is retained for the model
# comparisons.
#
# This script producees 500 .rds files in the autistic/ directory. These
# files are processed and aggregated into a final solution in the next
# script: 11-evaluate-best-growth-models_2024_01_05.R.

ncores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(ncores)
parallel::clusterSetRNGStream(cl)
invisible(parallel::clusterEvalQ(cl, {
    library(dplyr)
    library(netgrowr)
    source("./R/utils.R")
}))


labels <- map2(network_gv1, network_gv2, c)
formulas <- map2(f1, f2g, c) |>
    map2(labels, ~ {
        names(.x) <- .y
        return(.x)
    })

for (repetition in 1:500) {
    M <- map2(formulas, modelvars_rn, ~ {
        parallel::parLapply(
            cl = cl,
            X = .x,
            fun = netgrowr::mle_network_growth,
            data = .y,
            split_by = "month",
            label_with = "num_item_id"
        )
    })

    M <- map2(formulas, modelvars_rn, ~ {
        map(
            .x,
            netgrowr::mle_network_growth,
            data = .y,
            split_by = "month",
            label_with = "num_item_id",
            .progress = TRUE
        )
    })
    M1 <- map2(M, network_gv1, ~ {
        .x[.y]
    })

    M2g <- map2(M, network_gv2, ~ {
        .x[.y]
    })

    # Comparison 1: Each growth model (and combinations) over baseline ----
    X <- list()
    X[[1]] <- t(cbind(
        vapply(M1$autistic, netgrowr::model_comparison, numeric(7), restricted = M0$autistic),
        vapply(M1$nonautistic, netgrowr::model_comparison, numeric(7), restricted = M0$nonautistic)
    ))

    # Comparison 2: Preferential Acq. vs. LOA ----
    comparisons <- as.data.frame(rbind(
        c(full = network_gv2$autistic[3], restricted = network_gv1$autistic[2]),
        c(full = network_gv2$autistic[3], restricted = network_gv1$autistic[3]),
        c(full = network_gv2$autistic[2], restricted = network_gv1$autistic[1]),
        c(full = network_gv2$autistic[2], restricted = network_gv1$autistic[3]),
        c(full = network_gv2$autistic[1], restricted = network_gv1$autistic[2]),
        c(full = network_gv2$autistic[1], restricted = network_gv1$autistic[1]),
        c(full = network_gv2$nonautistic[3], restricted = network_gv1$nonautistic[2]),
        c(full = network_gv2$nonautistic[3], restricted = network_gv1$nonautistic[3]),
        c(full = network_gv2$nonautistic[2], restricted = network_gv1$nonautistic[1]),
        c(full = network_gv2$nonautistic[2], restricted = network_gv1$nonautistic[3]),
        c(full = network_gv2$nonautistic[1], restricted = network_gv1$nonautistic[2]),
        c(full = network_gv2$nonautistic[1], restricted = network_gv1$nonautistic[1])
    ))

    X[[2]] <- t(mapply(model_comp_helper, comparisons[["full"]], comparisons[["restricted"]], MoreArgs = list(M = c(M1$autistic, M1$nonautistic, M2g$autistic, M2g$nonautistic))))
    rownames(X[[2]]) <- apply(comparisons, 1, paste, collapse = "|")

    model_comparisons <- do.call("rbind", X)
    model_comparisons <- cbind(
        model_comparisons,
        p_fdr  = p.adjust(as.vector(model_comparisons[, "p"]), method = "fdr"),
        p_bonf = p.adjust(as.vector(model_comparisons[, "p"]), method = "bonferroni"),
        p_holm = p.adjust(as.vector(model_comparisons[, "p"]), method = "holm")
    )

    saveRDS(model_comparisons, file = sprintf("./autistic/model-comparisons-autistic_wCHILDES_%03d.rds", repetition - 1))
}
parallel::stopCluster(cl)
