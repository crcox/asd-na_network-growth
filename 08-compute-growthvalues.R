library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(igraph)
library(netgrowr)

# Load id key ----
d <- readRDS("./data/item_level_differences_bs_ci_bonf.rds")
m <- readRDS("./data/cdi-metadata.rds")
g <- upgrade_graph(readRDS("./network/child_net_graph.rds"))
assocnet <- as_adjacency_matrix(g, sparse = FALSE)

# N.B. all.equal(rowSums(assocnet), degree(g, mode = "out")) == TRUE
#     Thus, a row in `assocnet` represents the connection "from" the word
#     associated with the row "to" the word associated with each column.

vertex_ids <- tibble(
    vid = seq_len(vcount(g)),
    word = names(V(g))
)

vsoa_df <- d %>%
    left_join(vertex_ids, by = "word") %>%
    filter(!is.na(vid)) %>%
    select(vid, num_item_id, word, vsoa_nonautistic = na, vsoa_autistic = asd) %>%
    arrange(vid) %>%
    pivot_longer(
        cols = starts_with("vsoa_"),
        names_to = c("group"),
        names_prefix = "vsoa_",
        values_to = "vsoa"
    )

vsoa_df <- vsoa_df %>%
    mutate(
        by_20 = cut(vsoa, c(-Inf, seq(20, 660, by = 20), Inf), ordered_result = TRUE),
        by_40 = cut(vsoa, c(-Inf, seq(20, 660, by = 40), Inf), ordered_result = TRUE),
        by_60 = cut(vsoa, c(-Inf, seq(20, 660, by = 60), Inf), ordered_result = TRUE)
    )

vsoa_plot <- vsoa_df %>%
    mutate(across(starts_with("by_"), as.character)) %>%
    pivot_longer(
        cols = starts_with("by_"),
        names_to = "binsize",
        values_to = "bin"
    )

ggplot(vsoa_plot, aes(x = bin, fill = group)) +
    geom_bar(position = position_dodge()) +
    facet_wrap(~binsize, scales = "free_x")

vsoa_keys <- vsoa_df %>%
    group_by(group) %>%
    group_keys()

vsoa_lst <- vsoa_df %>%
    group_by(group) %>%
    group_split()

names(vsoa_lst) <- vsoa_keys$group

growthvalues <- map(vsoa_lst, function (vsoa, assocnet) {
    netgrowr::growth_values(
        assocnet,
        aoa_tbl = vsoa %>% select(word, by_20) %>% mutate(by_20 = as.numeric(by_20)),
        growth_models = c("preferential_attachment", "lure_of_the_associates", "preferential_acquisition")
    ) %>%
        left_join(
            vsoa %>% select(vocab_step = by_20) %>% mutate(month = as.numeric(vocab_step)) %>% distinct(),
            by = "month"
        ) %>%
        left_join(
            vsoa %>% select(word, vsoa, vsoa_bin = by_20) %>% distinct(),
            by = "word"
        ) %>%
        as_tibble()
}, assocnet = assocnet)

# Save growth values ----
saveRDS(growthvalues$autistic, file = "./network/growthvalues-autistic.rds")
saveRDS(growthvalues$nonautistic, file = "./network/growthvalues-nonautistic.rds")

# Save wide-form data for lexical growth modeling ----
# Specify and incorporate phonological baseline variables
phono_baseline <- m %>%
    select(
        word = cue_CoxHae,
        num_item_id,
        nphon,
        CHILDES_Freq,
        BiphonProb.avg,
        PNDC.avg
    ) %>%
    left_join(vertex_ids, by = "word") %>%
    filter(!is.na(vid))

growthvalues$autistic$group <- "autistic"
growthvalues$nonautistic$group <- "nonautistic"

modelvars <- growthvalues %>%
    list_rbind() %>%
    tidyr::pivot_wider(
        id_cols = c('word', 'month', 'learned', "vocab_step"),
        names_from = c(model, group),
        values_from = 'value'
    ) %>%
    left_join(phono_baseline, by = "word") %>%
    left_join(vsoa_df %>%
                  filter(group == "autistic") %>%
                  select(word, vsoa_autistic = vsoa, vsoa_bin_autistic = by_20),
              by = "word"
    ) %>%
    left_join(vsoa_df %>%
                  filter(group == "nonautistic") %>%
                  select(word, vsoa_nonautistic = vsoa, vsoa_bin_nonautistic = by_20),
              by = "word"
    )


saveRDS(modelvars, file = "./network/modelvars.rds")
