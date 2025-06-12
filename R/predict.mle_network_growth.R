predict.mle_network_growth <- function(model, dnew = NULL) {
    dnew <- if (is.null(dnew)) model$model
    netgrowr:::probability_node_added(
        coef(model),
        model$formula,
        dnew,
        split_by = "vocab_step"
    )
}
