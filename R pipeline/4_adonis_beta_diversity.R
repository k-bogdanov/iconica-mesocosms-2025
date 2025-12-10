# 16S
distance <- phyloseq::distance(exp1_ps_16s_r, "bray")
ado_df <- as(sample_data(exp1_ps_16s_r), "data.frame") %>%
  dplyr::select(P.level, Day)

# ado_df$Texture <- factor(ado_df$Texture, levels = c("Silt", "Silt Loam", "Sandy Loam"))

factors <- c(colnames(ado_df))

# create new vectors each time before the for loop
p_values <- vector(mode = "double", length = 0)
r_values <- vector(mode = "double", length = 0)

for(i in factors) {
  ado_res <- adonis2(distance ~ ado_df[[i]], ado_df)
  p_values = append(p_values, ado_res$`Pr(>F)`[1])
  r_values = append(r_values, ado_res$R2[1])
}

adonis_res_16s <- data.frame(factors, p_values, r_values) %>%
  mutate(Type = "16S") %>%
  dplyr::rename(Factor = factors, p = p_values, R2 = r_values)
adonis_res_16s 

# ITS
distance <- phyloseq::distance(exp1_ps_ITS_r, "bray")
ado_df <- as(sample_data(exp1_ps_ITS_r), "data.frame") %>%
  dplyr::select(P.level, Day)

factors <- c(colnames(ado_df))

# create new vectors each time before the for loop
p_values <- vector(mode = "double", length = 0)
r_values <- vector(mode = "double", length = 0)

for(i in factors) {
  ado_res <- adonis2(distance ~ ado_df[[i]], ado_df)
  p_values = append(p_values, ado_res$`Pr(>F)`[1])
  r_values = append(r_values, ado_res$R2[1])
}

adonis_res_its <- data.frame(factors, p_values, r_values) %>%
  mutate(Type = "ITS") %>%
  dplyr::rename(Factor = factors, p = p_values, R2 = r_values)
adonis_res_its

# apply(truesoil_all_metadata[4:21], FUN = shapiro.test, MARGIN = 2)

adonis_all <- bind_rows(adonis_res_16s, adonis_res_its) %>%
  arrange(Type, p, desc(R2))
adonis_all

# ordinate data using Non-metric multidimensional scaling (NMDS) on Bray–Curtis dissimilarity (distances)

set.seed(1509)
# 16S
exp1_ps_16s_ord <- ordinate(exp1_ps_16s_r, "NMDS", "bray")
exp1_ps_16s_ord

stressplot(exp1_ps_16s_ord)

stress_val_16s = grobTree(textGrob("SSP: p ≤ 0.001, R² = 0.139\nTime: p = 0.002, R² = 0.091", x = 0.50,  y = 0.10, hjust = 0, gp = gpar(col = "Black", fontsize = 8)))

exp1_16s_nmds <- plot_ordination(exp1_ps_16s_r, exp1_ps_16s_ord, justDF = TRUE) %>% 
  ggplot(aes(NMDS1, NMDS2)) +
  scale_color_manual(name = "SP (kg/ha)", values = c("#F3A26D", "#708A58")) +
  geom_point(aes(alpha = Day, color = P.level), size = 4) +
  my_theme + theme(legend.position = "none") +
  annotation_custom(stress_val_16s) +
  labs(color = "SP (kg/ha)", alpha = "Time (d)") +
  ggtitle("B")




# ITS
exp1_ps_ITS_ord <- ordinate(exp1_ps_ITS_r, "NMDS", "bray")
exp1_ps_ITS_ord

stressplot(exp1_ps_ITS_ord)

stress_val_ITS = grobTree(textGrob("SSP: p ≤ 0.001, R² = 0.224\nTime: p ≤ 0.001, R² = 0.153", x = 0.50,  y = 0.10, hjust = 0, gp = gpar(col = "Black", fontsize = 8)))

exp1_ITS_nmds <- plot_ordination(exp1_ps_ITS_r, exp1_ps_ITS_ord, justDF = TRUE) %>% 
  ggplot(aes(NMDS1, NMDS2)) +
  scale_color_manual(name = "SP (kg/ha)", values = c("#F3A26D", "#708A58")) +
  geom_point(aes(alpha = Day, color = P.level), size = 4) +
  my_theme + theme(legend.position = "right") +
  annotation_custom(stress_val_ITS) +
  labs(color = "SP (kg/ha)", alpha = "Time (d)") +
  ggtitle("D") +
  scale_x_continuous(limits = c(-1, 1)) +
  guides(color = guide_legend(title = expression("SSP (kg "*ha^{-1}*' '*year^{-1}*")"), override.aes = aes(label = ""))) +
  ylim(-1, 1)

legend <- get_legend(exp1_ITS_nmds)

exp1_ITS_nmds <- exp1_ITS_nmds + rremove("legend")
exp1_ITS_nmds