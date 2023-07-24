library(apaTables)
library(ez)

ttest_vec <- function(x) {
  c(
    tval = as.vector(x[['statistic']]),
    df = as.vector(x[['parameter']]),
    p = as.vector(x[['p.value']]),
    ci.l = as.vector(x[['conf.int']][1]),
    ci.u = as.vector(x[['conf.int']][2]),
    stderr = as.vector(x[['stderr']])
  )
}


load_to_list <- function(...) {
    files <- rlang::list2(...)
    purrr::map(files, function(filename) {
        e <- new.env()
        nm <- load(filename, envir = e)
        e[[nm]]
    })
}

# Load processed word associations ----
m <- readRDS("./data/cdi-metadata.rds")
growth_values <- list(
  autistic = readRDS("./network/growthvalues-autistic.rds"),
  nonautistic = readRDS("./network/growthvalues-nonautistic.rds")
)

# Combine into one dataframe ----
d <- growth_values %>%
    imap(~{
        .x$group <- factor(.y, levels = c("autistic", "nonautistic"))
        return(.x)
    }) %>%
    list_rbind() %>%
    filter(learned == TRUE) %>%
    droplevels()

# Report means and one-way t-tests against zero ----
means_tbl <- tapply(d$zscore, d[c("model", "group")], mean)
print(means_tbl, digits = 3)

stdev_tbl <- tapply(d$zscore, d[c("model", "group")], sd)
print(means_tbl, digits = 3)

barplot(t(means_tbl[c(1,3,2), ]), beside = TRUE)

ttests_list <- lapply(split(d$zscore, list(d$model, d$group)), t.test)
ttests_tbl <- t(vapply(ttests_list, ttest_vec, numeric(6)))
print(ttests_tbl, digits = 3)


shared_words <- intersect(unique(d$word[d$group=="autistic"]), unique(d$word[d$group=="nonautistic"]))

d_aov <- d %>%
    filter(word %in% shared_words) %>%
    mutate(word = as.factor(word))

eza <- ezANOVA(data = d_aov, wid = .(word), dv = .(zscore), within = .(model, group))
apa.ezANOVA.table(eza)

# Simple effects by model ----
## Paired t-test
simple_effects_by_model <- lapply(split(d_aov, d_aov$model), function(x) {t.test(zscore ~ group, data = x, paired = TRUE)})
tt_by_model <- t(vapply(simple_effects_by_model, ttest_vec, numeric(6)))
print(tt_by_model, 3)

w <- pivot_wider(d_aov, id_cols = c("model", "word", "month", "aoa"), names_from = "group", values_from = "zscore")
lapply(split(w, w$model), function(x) abs(mean(x$autistic - x$nonautistic, na.rm = TRUE) / sd(x$autistic - x$nonautistic, na.rm = TRUE)))

# Simple effects by condition ----
## Repeated measures ANOVA
apaTables::apa.ezANOVA.table(
  ez::ezANOVA(subset(d_aov, group == "autistic"),
              dv = .(zscore),
              wid = .(word),
              within = .(model)
  )
)

apaTables::apa.ezANOVA.table(
  ez::ezANOVA(subset(d_aov, group == "nonautistic"),
              dv = .(zscore),
              wid = .(word),
              within = .(model)
  )
)

## Simple effects diving down another level.
t.test(zscore ~ model, data = droplevels(subset(d_aov, model %in% c("preferential_acquisition", "lure_of_the_associates") & group == "autistic")), paired = TRUE)
t.test(zscore ~ model, data = droplevels(subset(d_aov, model %in% c("preferential_acquisition", "lure_of_the_associates") & group == "nonautistic")), paired = TRUE)
