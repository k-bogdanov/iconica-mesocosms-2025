# 16S

alpha_sum_16s <- estimate_richness(exp1_ps_16s_r, measures = c("Observed", "Shannon"))

# this calculates eveness and adds it to the prior file

Evenness_16s <- microbiome::evenness(exp1_ps_16s_r, "pielou")
alpha_sum_16s$Pielou <- Evenness_16s$pielou

# combine results with sample_data. This is useful if we wish to plot later against other variables or information.
alpha_meta_16s <- data.frame(alpha_sum_16s, sample_data(exp1_ps_16s_r))

# convert data from wide to long for easier plotting
alpha_meta_16s_long = gather(alpha_meta_16s, Alpha, alpha_value, c("Observed", "Shannon", "Pielou"), factor_key = FALSE)

observed_sum_16s <- alpha_meta_16s_long %>%
  group_by(Alpha, P.level, Day) %>%
  dplyr::summarise(
    N = dplyr::n(),
    median = median(alpha_value, na.rm = TRUE), # replace with mean
    IQR = IQR(alpha_value, na.rm = TRUE), # replace with sd
    .groups = "drop"
  )


# making alpha diversity plot

alpha_plot_16s <- ggplot(observed_sum_16s,
                         aes(x = Day, y = median, colour = P.level, group = P.level)) +
  geom_point(size = 3) +
  geom_line(data = observed_sum_16s,
            aes(x = Day, y = median, group = P.level, color = P.level)) +
  geom_errorbar(aes(ymin = median - observed_sum_16s$IQR, # or mean - sd
                    ymax = median + observed_sum_16s$IQR), # or mean + sd
                width = 0.1) +
  my_theme +
  scale_color_manual(values = c("#F3A26D", "#708A58")) +
  scale_shape_manual(values = c(4, 18, 19)) +
  labs(x = "Time (d)", y = "Alpha diversity", color = "SP (kg/ha)") +
  facet_grid(Alpha~., scales = "free_y") +
  theme(legend.position = "right",
        axis.ticks = element_blank()) +
  ggtitle("A")

# check normality
columns <- c("Observed", "Shannon", "Pielou")

for (col in columns) {
  print(shapiro_test(alpha_meta_16s[[col]]))
  print(ggdensity(alpha_meta_16s[[col]], title = paste("Density Plot -", col)))
}

# stats by P level
# wilcox test
columns <- c("Observed", "Shannon", "Pielou")
wilcox_test_results <- list()

for (col in columns) {
  result <- alpha_meta_16s %>%
    wilcox_test(as.formula(paste(col, "~ P.level"))) %>%
    mutate(variable = col)
  
  wilcox_test_results[[col]] <- result
}

bind_rows(wilcox_test_results)


# stats by day
for (col in columns) {
  result <- alpha_meta_16s %>%
    group_by(P.level) %>%
    t_test(as.formula(paste(col, "~ Day"))) %>%
    filter(p.adj < 0.05) %>%
    select(-statistic, -n1, -n2, -statistic)
  
  t_test_results[[col]] <- result
}


bind_rows(t_test_results)

# ITS

alpha_sum_ITS <- estimate_richness(exp1_ps_ITS_r, measures = c("Observed", "Shannon"))

# this calculates eveness and adds it to the prior file

Evenness_ITS <- microbiome::evenness(exp1_ps_ITS_r, "pielou")
alpha_sum_ITS$Pielou <- Evenness_ITS$pielou

# combine results with sample_data. This is useful if we wish to plot later against other variables or information.
alpha_meta_ITS <- data.frame(alpha_sum_ITS, sample_data(exp1_ps_ITS_r))

# convert data from wide to long for easier plotting
alpha_meta_ITS_long = gather(alpha_meta_ITS, Alpha, alpha_value, c("Observed", "Shannon", "Pielou"), factor_key = FALSE)

# summarise data
observed_sum_ITS <- alpha_meta_ITS_long %>%
  group_by(Alpha, P.level, Day) %>%
  dplyr::summarise(
    N = dplyr::n(),
    median = median(alpha_value, na.rm = TRUE),
    IQR = IQR(alpha_value, na.rm = TRUE),
    mean = mean(alpha_value, na.rm = TRUE),
    std = std(alpha_value),
    .groups = "drop"
  )


# making alpha diversity plot

alpha_plot_ITS <- ggplot(observed_sum_ITS,
                         aes(x = Day, y = mean, colour = P.level, group = P.level)) +
  geom_point(size = 3) +
  geom_line(data = observed_sum_ITS,
            aes(x = Day, y = mean, group = P.level, color = P.level)) +
  geom_errorbar(aes(ymin = mean - observed_sum_ITS$std, 
                    ymax = mean + observed_sum_ITS$std), 
                width = 0.1) +
  my_theme +
  scale_color_manual(values = c("#F3A26D", "#708A58")) +
  labs(x = "Time (d)", y = "Alpha diversity", color = "SP (kg/ha)") +
  facet_grid(Alpha~., scales = "free_y") +
  theme(legend.position = "right",
        axis.ticks = element_blank()) +
  ggtitle("C")
alpha_plot_ITS

# legend <- get_legend(alpha_plot_ITS)

alpha_plot_ITS_exp1 <- alpha_plot_ITS + rremove("legend")
alpha_plot_16s_exp1 <- alpha_plot_16s + rremove("legend")


# check normality
columns <- c("Observed", "Shannon", "Pielou")

for (col in columns) {
  print(shapiro_test(alpha_meta_ITS[[col]]))
  print(ggdensity(alpha_meta_ITS[[col]], title = paste("Density Plot -", col)))
}

# stats by P level
# wilcox test
columns <- c("Observed")
wilcox_test_results <- list()

for (col in columns) {
  result <- alpha_meta_ITS %>%
    wilcox_test(as.formula(paste(col, "~ P.level"))) %>%
    mutate(variable = col)
  
  wilcox_test_results[[col]] <- result
}

bind_rows(wilcox_test_results)



# t test
columns <- c("Shannon", "Pielou")
t_test_results <- list()

for (col in columns) {
  result <- alpha_meta_ITS %>%
    group_by(Day) %>%
    t_test(as.formula(paste(col, "~ P.level"))) %>%
    mutate(variable = col) %>%
    filter(p < 0.05)
  
  t_test_results[[col]] <- result
}

bind_rows(t_test_results)


# stats by day
for (col in columns) {
  result <- alpha_meta_ITS %>%
    group_by(P.level) %>%
    t_test(as.formula(paste(col, "~ Day"))) %>%
    filter(p < 0.05) %>%
    select(-statistic, -n1, -n2, -statistic)
  
  t_test_results[[col]] <- result
}


bind_rows(t_test_results)