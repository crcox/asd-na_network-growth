library(dplyr)
library(purrr)
library(netgrowr)
library(parallel)
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
modelvars_df <- readRDS("./network/modelvars_vsoa_2025_05_14.rds") %>%
    mutate(label = paste(group, num_item_id))

modelvars <- modelvars_df %>%
    group_by(group) %>%
    group_split()

names(modelvars) <- modelvars_df %>%
    group_by(group) %>%
    group_keys() %>%
    pull(group)

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


# Define models by constructing formulas ----

# Empty (random selection) model
# +++ All words assigned equal probability
fE <- formula(aoa ~ 1)
ME <- netgrowr::mle_network_growth(fE, data = na.omit(modelvars_df), split_by = "month", label_with = "label")

# Psycholinguistic baseline model
# +++ number of phonemes
# +++ CHILDES Frequency
# +++ phonotactic probability (biphone)
# +++ phonological neighborhood density (KU child corpus)
f0 <- update(fE, ~ . + I(Z(log(nphon + 1))) + I(Z(log(CHILDES_Freq + 1))) + Z(BiphonProb.avg) + Z(PNDC.avg))
M0 <- map(modelvars, ~ {
    netgrowr::mle_network_growth(f0, data = na.omit(.x), split_by = "month", label_with = "num_item_id")
})

# M0gR <- netgrowr::mle_network_growth(f0, data = na.omit(modelvars_df), split_by = "month", label_with = "label")
# f0g <- update(f0, ~ (.) * group)
# M0gF <- netgrowr::mle_network_growth(f0g, data = na.omit(modelvars_df), split_by = "month", label_with = "label")

# Networks over control
network_gv1 <- list(
    autistic = c(
        "preferential_attachment_child_autistic + preferential_attachment_childes_autistic",
        "preferential_acquisition_child_autistic + preferential_acquisition_childes_autistic",
        "lure_of_the_associates_child_autistic + lure_of_the_associates_childes_autistic"
    ),
    nonautistic = c(
        "preferential_attachment_child_nonautistic + preferential_attachment_childes_nonautistic",
        "preferential_acquisition_child_nonautistic + preferential_acquisition_childes_nonautistic",
        "lure_of_the_associates_child_nonautistic + lure_of_the_associates_childes_nonautistic"
    )
)
f1 <- map_depth(network_gv1, 2, append_formula, f = f0, .is_node = list_or_vec)


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
formulas <- map2(f1, f2g, c) %>%
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

    saveRDS(model_comparisons, file = sprintf("./autistic/20250514/model-comparisons-autistic_wCHILDES_%03d.rds", repetition - 1))
}
parallel::stopCluster(cl)
