# 16S

exp1_ps_16s_df <- exp1_ps_16s_r %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()

# stats
wilcox_16s_P_level <- exp1_ps_16s_df %>%
  group_by(Species) %>%
  filter(n_distinct(P.level) > 1) %>%
  wilcox_test(Abundance ~ P.level, paired = FALSE) %>%
  adjust_pvalue(method = "fdr") %>% 
  filter(p.adj < 0.05)
wilcox_16s_P_level

# combine significantly affected taxa with initial table
only_sig_16S <- exp1_ps_16s_df %>% 
  filter(Species %in% wilcox_16s_P_level$Species) %>%
  filter(!str_detect(Genus, "Unidentified")) %>%
  mutate(across(Species, ~ ifelse(grepl("Unidentified", Species), "sp", .))) %>%
  mutate(across(Genus, ~ ifelse(grepl("Burkholderia-Caballeronia-Paraburkholderia", Genus), "Paraburkholderia", .))) %>%
  filter(!Genus %in% c("Clostridium", "Terrimonas", "Pedobacter", "Pelotalea"))

# plot
sig_affected_taxa_16S <- ggplot(only_sig_16S, aes(
  x = Day,
  y = Abundance
)) +
  geom_boxplot(aes(fill = P.level)) +
  facet_wrap(Genus ~ Species, scales = "free_y") +
  scale_y_log10(labels = scales::scientific) +
  scale_fill_manual(name = expression("SSP (kg "*ha^{-1}*")"), values = c("#F3A26D", "#708A58")) +
  labs(x = "Time (d)", y = "Relative abundance (log10)") +
  theme(axis.title = element_text(face = "bold", size = 16),
        panel.background = element_rect(fill = "white", colour = "grey90"),
        strip.text = element_text(face = "italic", size = 12),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1, "cm"),
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_blank(),
        panel.grid = element_line(color = "grey90"),
        panel.grid.major = element_blank(),
        plot.margin = margin(l = 10, b = 10, r = 10, t = 10),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(colour = "grey90", fill ="white"),
        legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key = element_rect(fill = NA))

# ITS

exp1_ps_ITS_df <- exp1_ps_ITS_r %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()

# stats
wilcox_ITS_P_level <- exp1_ps_ITS_df %>%
  group_by(Species) %>%
  filter(n_distinct(P.level) > 1) %>%
  wilcox_test(Abundance ~ P.level, paired = FALSE) %>%
  adjust_pvalue(method = "fdr") %>% 
  filter(p.adj < 0.05)
wilcox_ITS_P_level

# combine significantly affected taxa with initial table
only_sig_ITS <- exp1_ps_ITS_df %>% 
  filter(Species %in% wilcox_ITS_P_level$Species) %>%
  mutate(across(Species, ~ ifelse(grepl("Unidentified", Species), "sp", .))) %>%
  mutate(across(Genus, ~ifelse(grepl("_gen_Incertae_sedis", Genus), "Unidentified", .))) %>%
  filter(!str_detect(Genus, "Unidentified"), !Genus %in% c("Linnemannia"))

# plot
sig_affected_taxa_ITS <- ggplot(only_sig_ITS, aes(
  x = Day,
  y = Abundance
)) +
  geom_boxplot(aes(fill = P.level)) +
  facet_wrap(Genus ~ Species, scales = "free_y") +
  scale_y_log10(labels = scales::scientific) +
  scale_fill_manual(name = expression("SSP (kg "*ha^{-1}*")"), values = c("#F3A26D", "#708A58")) +
  labs(x = "Time (d)", y = "Relative abundance (log10)") +
  theme(axis.title = element_text(face = "bold", size = 16),
        panel.background = element_rect(fill = "white", colour = "grey90"),
        strip.text = element_text(face = "italic", size = 12),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1, "cm"),
        axis.text = element_text(color = "grey40"),
        axis.ticks = element_blank(),
        panel.grid = element_line(color = "grey90"),
        panel.grid.major = element_blank(),
        plot.margin = margin(l = 10, b = 10, r = 10, t = 10),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        strip.background = element_rect(colour = "grey90", fill ="white"),
        legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key = element_rect(fill = NA))