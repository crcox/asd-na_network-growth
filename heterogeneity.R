library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

vsoa <- readRDS("data/vsoa-autistic-nonautistic-ndar-id-fix-remodel-v2.rds")

models <- map(1:680, ~{
    readRDS("")
})

cdi <- readRDS("../asd-na_item-bs/data/asd_na-osg-2025-05-20.rds")

d <- left_join(
    cdi,
    vsoa |> select(group, num_item_id, vsoa)
)

q <- d |>
    group_by(group, subjectkey) |>
    summarize(
        prop = mean(produced == (vsoa <= nproduced)),
        tpr  = mean( produced &  (vsoa <= nproduced)),
        fpr  = mean(!produced &  (vsoa <= nproduced)),
        tnr  = mean(!produced & !(vsoa <= nproduced)),
        fnr  = mean( produced & !(vsoa <= nproduced))
    ) |>
    mutate(
        adjacc = (tpr + tnr) / 2,
        dprime = qnorm(tpr) - qnorm(fpr)
    )

q |>
    group_by(group) |>
    summarize(
        y = mean(adjacc),
        se = sd(adjacc) / sqrt(n()),
        ci = se * qt(0.025, n() - 1),
        ymin = y - ci,
        ymax = y + ci
    ) |>
    ggplot(aes(x = group, y = y)) +
    geom_pointrange(aes(ymin = ymin, ymax = ymax))

t.test(adjacc ~ group, data = q)
