library(dplyr)
library(purrr)
library(tidyr)
library(netgrowr)
source("./R/utils.R")


# Note: 15 May 2025
# Including a group interaction might help determine if growth models (and
# confounding variables) are weighted significantlt differently between groups.


model_comp_helper <- function(full, restricted, M) {
    return(netgrowr::model_comparison(M[[full]], M[[restricted]]))
}


# Load data ----
# Month 30 is the highest month we have AoA data for.
# If we model month 30, all unknown words will be learned.
# To avoid this, month 30 is not modeled.
modelvars <- readRDS("./network/modelvars-vsoa-20250514.rds")
modelvars <- modelvars %>%
    drop_na() %>%
    mutate(
        label = paste(
            if_else(group=="autistic", "a", "n"),
            month,
            num_item_id,
            sep = "_"
        ),
        group = factor(group, levels = c("autistic", "nonautistic"))
    )


# Define models by constructing formulas ----

# Empty (random selection) model
# +++ All words assigned equal probability
fE <- formula(aoa ~ 1)
ME <- map(c(autistic = "autistic", nonautistic = "nonautistic"), function(f, g, .data) {
    netgrowr::mle_network_growth(
        f,
        data = .data %>% drop_na() %>% filter(group == g, model == "preferential_acquisition"),
        split_by = "month",
        label_with = "label"
    )
}, f = fE, .data = modelvars)
names(ME) <- c("autistic", "nonautistic")
map2(M0, ME, ~{
    model_comparison(.x, .y)
}) %>%
    map(~{as_tibble(as.list(.))}) %>%
    list_rbind() %>%
    bind_cols(tibble(group = c("autistic", "nonautistic")))

# Psycholinguistic baseline model
# +++ number of phonemes
# +++ CHILDES Frequency
# +++ phonotactic probability (biphone)
# +++ phonological neighborhood density (KU child corpus)
f0 <- update(fE, ~ . + I(Z(log(nphon + 1))) + I(Z(log(CHILDES_Freq + 1))) + Z(BiphonProb.avg) + Z(PNDC.avg))
M0 <- map(c(autistic = "autistic", nonautistic = "nonautistic"), function(f, g, .data) {
    netgrowr::mle_network_growth(
        f,
        data = .data %>% drop_na() %>% filter(group == g, model == "preferential_acquisition"),
        split_by = "month",
        label_with = "label"
    )
}, f = f0, .data = modelvars)
names(M0) <- c("autistic", "nonautistic")

x <- modelvars %>%
    select(num_item_id, nphon, CHILDES_Freq, BiphonProb.avg, PNDC.avg) %>%
    distinct()

q <- model.matrix(update(f0, ~ . -1), mutate(x, aoa = 1))
r <- cor(q)
rr <- cor(t(q))
s <- svd(rr)
v <- varimax(s$u[,1:2])
#p <- psych::pca(rr, nfactors = 2)
x <- bind_cols(x, tibble(RC1 = v$loadings[,1], RC2 = v$loadings[,2]))

modelvars <- left_join(modelvars, select(x, num_item_id, starts_with("RC")))


f0RC <- update(fE, ~ . + RC1 + RC2)
M0RC <- map(c(autistic = "autistic", nonautistic = "nonautistic"), function(f, g, .data) {
    netgrowr::mle_network_growth(
        f,
        data = .data %>% drop_na() %>% filter(group == g, model == "preferential_acquisition"),
        split_by = "month",
        label_with = "label"
    )
}, f = f0RC, .data = modelvars)
names(M0RC) <- c("autistic", "nonautistic")
map2(M0RC, ME, ~{
    model_comparison(.x, .y)
}) %>%
    map(~{as_tibble(as.list(.))}) %>%
    list_rbind() %>%
    bind_cols(tibble(group = c("autistic", "nonautistic")))


# Networks over control
f1 <- update(f0, ~ . + value)
x <- expand_grid(
    model = c("preferential_acquisition", "lure_of_the_associates"),
    group = c("autistic", "nonautistic")
)
M1 <- map2(x$model, x$group, function(m, g, .formula, .data) {
    netgrowr::mle_network_growth(
        .formula,
        data = .data %>% drop_na() %>% filter(group == g, model == m),
        split_by = "month",
        label_with = "label"
    )
}, .formula = f1, .data = modelvars)
names(M1) <- paste(x$model, x$group, sep = '_')

list(
    model_comparison(M1$preferential_acquisition_autistic, M0$autistic),
    model_comparison(M1$preferential_acquisition_nonautistic, M0$nonautistic),
    model_comparison(M1$lure_of_the_associates_autistic, M0$autistic),
    model_comparison(M1$lure_of_the_associates_nonautistic, M0$nonautistic)
) %>% map(~{as_tibble(as.list(.))}) %>% list_rbind() %>% bind_cols(x)


f1RC <- update(f0RC, ~ . + value)
x <- expand_grid(
    model = c("preferential_acquisition", "lure_of_the_associates"),
    group = c("autistic", "nonautistic")
)
M1RC <- map2(x$model, x$group, function(m, g, .formula, .data) {
    netgrowr::mle_network_growth(
        .formula,
        data = .data %>% drop_na() %>% filter(group == g, model == m),
        split_by = "month",
        label_with = "label"
    )
}, .formula = f1RC, .data = modelvars)
names(M1RC) <- paste(x$model, x$group, sep = '_')

list(
    model_comparison(M1RC$preferential_acquisition_autistic, M0RC$autistic),
    model_comparison(M1RC$preferential_acquisition_nonautistic, M0RC$nonautistic),
    model_comparison(M1RC$lure_of_the_associates_autistic, M0RC$autistic),
    model_comparison(M1RC$lure_of_the_associates_nonautistic, M0RC$nonautistic)
) %>% map(~{as_tibble(as.list(.))}) %>% list_rbind() %>% bind_cols(x)



fInt <- update(f1, ~ (.) * group)
MInt <- map(
    c(
        preferential_acquisition = "preferential_acquisition",
        lure_of_the_associates = "lure_of_the_associates"
    ),
    function(m, .formula, .data) {
        netgrowr::mle_network_growth(
            .formula,
            data = .data %>% drop_na() %>% filter(model == m),
            split_by = "month",
            label_with = "label"
        )
    }, .formula = fInt, .data = modelvars
)


fGroupRC <- update(f1RC, ~ . + group)
MGroupRC <- map(
    c(
        preferential_acquisition = "preferential_acquisition",
        lure_of_the_associates = "lure_of_the_associates"
    ),
    function(m, .formula, .data) {
        netgrowr::mle_network_growth(
            .formula,
            data = .data %>% drop_na() %>% filter(model == m),
            split_by = "month",
            label_with = "label"
        )
    }, .formula = fGroupRC, .data = modelvars
)
tmp <- solve(MGroupRC$preferential_acquisition$hessian)

fIntRC <- update(f1RC, ~ (.) * group)
MIntRC <- map(
    c(
        preferential_acquisition = "preferential_acquisition",
        lure_of_the_associates = "lure_of_the_associates"
    ),
    function(m, .formula, .data) {
        netgrowr::mle_network_growth(
            .formula,
            data = .data %>% drop_na() %>% filter(model == m),
            split_by = "month",
            label_with = "label"
        )
    }, .formula = fIntRC, .data = modelvars
)
tmp <- solve(MIntRC$preferential_acquisition$hessian)


map(MIntRC, ~ {
    se <- sqrt(diag(solve(.x$hessian[-1, -1])))
    x <- coef(.x)[-1]
    tval <- x / se
    return(tval)
}) %>%
    map(~{as_tibble(as.list(.))}) %>%
    list_rbind() %>%
    bind_cols(tibble(model = c("pref. acq.", "lure of assoc.")))

MCombRC <- map(
    c(
        preferential_acquisition = "preferential_acquisition",
        lure_of_the_associates = "lure_of_the_associates"
    ),
    function(m, .formula, .data) {
        netgrowr::mle_network_growth(
            .formula,
            data = .data %>% drop_na() %>% filter(model == m),
            split_by = "month",
            label_with = "label"
        )
    }, .formula = f1RC, .data = modelvars
)

list(
    model_comparison(MIntRC$preferential_acquisition, MCombRC$preferential_acquisition),
    model_comparison(MIntRC$lure_of_the_associates, MCombRC$lure_of_the_associates)
) %>% map(~{as_tibble(as.list(.))}) %>% list_rbind() %>% bind_cols(tibble(model = c("pref. acq.", "lure of assoc.")))
list(
    model_comparison(MGroupRC$preferential_acquisition, MCombRC$preferential_acquisition),
    model_comparison(MGroupRC$lure_of_the_associates, MCombRC$lure_of_the_associates)
) %>% map(~{as_tibble(as.list(.))}) %>% list_rbind() %>% bind_cols(tibble(model = c("pref. acq.", "lure of assoc.")))
# Define cluster for parallel computing and optimize models ----

# The optimization process is not guaranteed to find a global minimum.
# Different runs of the optimization will yield numeric differences, and
# rarely very bad solutions are obtained. Thus, the optimization process is
# repeated 500 times for every model fit. The minimum negative log likelihood
# for each model over the 500 optimization attempts is retained for the model
# comparisons.
#
# In each repetition, 22 models are optimized. The process can be greatly
# accellerated with parallel computing. We define the 22 formulas defining
# each model in a list, then use parallel::parLapply to fit the models in
# parallel.
#
# This script producees 500 Rdata files in the network/ directory. These
# files are processed and aggregated into a final solution in the next
# script: 05-aggregate-growth-models.R.

ncores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(ncores)
parallel::clusterSetRNGStream(cl)
invisible(parallel::clusterEvalQ(cl, {
    library('dplyr')
    library('netgrowr')
    source('./R/utils.R')
}))


for (repetition in 1:500) {
    M <- parallel::parLapply(
        cl = cl,
        X = c(f1, f2n, f3n, f2g),
        fun = netgrowr::mle_network_growth,
        data = modelvars,
        split_by = "month",
        label_with = "num_item_id"
    )
    names(M) <- c(
        network_gv1,
        network_gv2,
        network_gv3,
        autistic_gv2,
        nonautistic_gv2
    )
    m1 <- m[network_gv1]
    m2n <- m[network_gv2]
    m3n <- m[network_gv3]
    m2g <- m[c(autistic_gv2, nonautistic_gv2)]

    # comparison 1: each growth model (and combinations) over baseline ----
    x <- list()
    x[[1]] <- t(vapply(
            c(m1, m2n, m3n),
            netgrowr::model_comparison,
            numeric(7),
            restricted = m0
    ))

    # comparison 2: preferential acq. vs. loa ----
    comparisons <- as.data.frame(rbind(
        c(full = 'preferential_acquisition_autistic+lure_of_the_associates_autistic', restricted = 'preferential_acquisition_autistic'),
        c(full = 'preferential_acquisition_autistic+lure_of_the_associates_autistic', restricted = 'lure_of_the_associates_autistic'),
        c(full = 'preferential_acquisition_nonautistic+lure_of_the_associates_nonautistic', restricted = 'preferential_acquisition_nonautistic'),
        c(full = 'preferential_acquisition_nonautistic+lure_of_the_associates_nonautistic', restricted = 'lure_of_the_associates_nonautistic')
    ))
    x[[2]] <- t(mapply(
        model_comp_helper,
        comparisons[["full"]],
        comparisons[["restricted"]],
        moreargs = list(m = c(m1, m2g))
    ))
    rownames(x[[2]]) <- apply(comparisons, 1, paste, collapse = "|")

    # comparison 3: autistic vs. nonautistic (preferential acq. only) ----
    comparisons <- as.data.frame(rbind(
        c(full = 'preferential_acquisition_autistic+preferential_acquisition_nonautistic', restricted = 'preferential_acquisition_autistic'),
        c(full = 'preferential_acquisition_autistic+preferential_acquisition_nonautistic', restricted = 'preferential_acquisition_nonautistic')
    ))
    x[[3]] <- t(mapply(
        model_comp_helper,
        comparisons[["full"]],
        comparisons[["restricted"]],
        moreargs = list(m = c(m1, m2n))
    ))
    rownames(x[[3]]) <- gsub(
        "preferential_acquisition_",
        "",
        apply(comparisons, 1, paste, collapse = "|")
    )

    model_comparisons <- do.call('rbind', X)
    model_comparisons <- cbind(
        model_comparisons,
        p_fdr  = p.adjust(as.vector(model_comparisons[, 'p']), method = 'fdr'),
        p_bonf = p.adjust(as.vector(model_comparisons[, 'p']), method = 'bonferroni'),
        p_holm = p.adjust(as.vector(model_comparisons[, 'p']), method = 'holm')
    )
    save(
        model_comparisons,
        file = sprintf('./network/model-comparisons-%03d.Rdata', repetition - 1)
    )
}
parallel::stopCluster(cl)
