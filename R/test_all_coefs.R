likelihood_ratio <- function(neg_loglik_full, neg_loglik_restricted) {
    return(2 * (neg_loglik_restricted - neg_loglik_full))
}

loglikelihood <- function(p) {
    if (is.matrix(p)) {
        return(colSums(log(p), na.rm = TRUE))
    } else {
        return(sum(log(p), na.rm = TRUE))
    }
}

ratio_of_strengths <- function (data, formula, beta, hold_const = NULL)
{
    learned <- data$learned
    unknown <- data$unknown
    d <- model.matrix(formula, data, na.action = "na.fail")
    if (!is.null(hold_const)) {
        z <- startsWith(colnames(d), hold_const)
        d[, z] <- mean(d[, z])
    }
    x <- exp(d %*% as.matrix(beta))
    return(t(t(x[learned, , drop = FALSE])/colSums(x[unknown,
                                                     , drop = FALSE])))
}

probability_node_added <- function (beta, formula, data, split_by, label_with = NULL, hold_const = NULL)
{
    d <- na.omit(netgrowr:::get_all_vars_with_split_and_labels(formula,
                                                    split_by, label_with, data))
    d$learned <- d[[1]] == d[[split_by]]
    d$unknown <- d[[1]] >= d[[split_by]]
    X <- split(d, d[[split_by]])
    if (!is.null(label_with)) {
        X <- lapply(X, function(x) {
            rownames(x) <- x[[label_with]]
            return(x)
        })
    }
    p <- do.call("rbind", lapply(X, ratio_of_strengths,
                                 beta = beta, formula = formula, hold_const = hold_const))
    if (ncol(p) == 1) {
        labs <- rownames(p)
        p <- as.vector(p)
        names(p) <- labs
    }
    return(p)
}

test_all_coefs <- function(m, split_by = "vocab_step") {
    f <- formula(vsoa_bin ~ (RC1 + RC2 + RC3 + loa1_z + pac1_z) * group)
    d <- model.frame(m)
    p <- probability_node_added(
        beta = coef(m),
        formula = f,
        data = d,
        split_by = split_by
    )
    x_full <- -loglikelihood(p)
    params <- names(coef(m))[-1]
    x_restricted <- map_dbl(params, ~{
        p <- probability_node_added(
            beta = coef(m),
            formula = f,
            data = d,
            split_by = split_by,
            hold_const = .x
        )
        -loglikelihood(p)
    })
    chisq <- map_dbl(x_restricted, ~ likelihood_ratio(x_full, .x))
    pval <- pchisq(chisq, df = 1, lower.tail = FALSE)
    BIC <- chisq - log(m$nobs)
    tibble(
        params,
        coef = coef(m)[params],
        df = 1,
        nobs = m$nobs,
        L0 = x_restricted,
        L1 = x_full,
        chisq = chisq,
        p = pval,
        BIC = BIC
    )
}

test_all_coefs(m)
