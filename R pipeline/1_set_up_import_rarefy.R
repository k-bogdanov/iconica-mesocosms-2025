# set up a theme
my_theme <- theme(axis.title.y = element_text(size = 12, face = "bold", vjust = 3),
                  axis.title.x = element_text(size = 12, face = "bold", vjust = -1),
                  panel.background = element_rect(fill = "white", colour = "grey90"),
                  strip.text.x = element_text(size = 12, colour = 'black'),
                  legend.key.width = unit(0.5, "cm"),
                  legend.key.height = unit(0.5, "cm"),
                  axis.text = element_text(color = "grey40"),
                  axis.ticks = element_blank(),
                  strip.text = element_text(colour = "black", size = 12),
                  panel.grid = element_line(color = "grey90"),
                  panel.grid.major = element_blank(),
                  plot.margin = margin(l = 5, b = 5, r = 5, t = 5),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold"),
                  strip.background = element_rect(colour = "grey90", fill ="white"),
                  legend.position = "right", legend.title = element_text(size = 12),
                  legend.text = element_text(size = 12),
                  legend.key = element_rect(fill = NA))

# import phyloseq objects
exp1_ps_16s <- readRDS("~/iconica_exp1_16s_phyloseq.rds") %>%
  name_na_taxa(na_label = "Unidentified <tax> (<rank>)")
exp1_ps_ITS <- readRDS("~/iconica_exp1_ITS_phyloseq.rds") %>%
  name_na_taxa(na_label = "Unidentified <tax> (<rank>)")


# remove any samples with less than 100 sequences
# exp1_ps_ITS_p <- prune_samples(sample_sums(exp1_ps_ITS) >= 100, exp1_ps_ITS)

# calculate rarefaction curves, then plot and color by subsite
# ggrare(exp1_ps_ITS_p, step = 1000, se = FALSE, parallel = TRUE) +
#  scale_x_continuous(limits = c(1, 65000)) +
#  geom_vline(xintercept = 5000)

set.seed(1509)

exp1_ps_16s_r <- rarefy_even_depth(exp1_ps_16s, 10000, rngseed = TRUE)
exp1_ps_ITS_r <- rarefy_even_depth(exp1_ps_ITS, 7000, rngseed = TRUE)


tax_table(exp1_ps_ITS_r)[, colnames(tax_table(exp1_ps_ITS_r))] <- gsub(tax_table(exp1_ps_ITS_r)[, colnames(tax_table(exp1_ps_ITS_r))], pattern = "[a-z]__", replacement = "")

# define factors
exp1_ps_16s_r <- exp1_ps_16s_r %>%
  ps_mutate(Day = factor(Day, levels = c("1", "2", "3", "4", "7")),
            P.level = factor(P.level, levels = c("0", "188")))

exp1_ps_ITS_r <- exp1_ps_ITS_r %>%
  ps_mutate(Day = factor(Day, levels = c("1", "2", "3", "4", "7")),
            P.level = factor(P.level, levels = c("0", "188")))

# export metadata
metadata_exp_1 <- data.frame(sample_data(exp1_ps_16s_r))
