# Setting up colors
p0 <- "#F3A26D"
p188 <- "#708A58"

library(common)
# Setting up labels for facets
Ntreatments <- as_labeller(c("14NH4 14NO3" = "NH[4]*''*NO[3]",
                             "15NH4 14NO3" = "''^15*NH[4]*''*NO[3]",
                             "14NH4 15NO3" = "''*NH[4]*''^15*NO[3]",
                             "14NH414NO3" = "''*NH[4]*''*NO[3]",
                             "15NH414NO3" = "''^15*NH[4]*''*NO[3]",
                             "14NH415NO3" = "''*NH[4]*''^15*NO[3]"), default = label_parsed)

OMtreatments <- as_labeller(c("15N-OM" = "''^13*C*' '*'+'*''^15*N*' '*roots",
                              "14N-OM" = "Control*' '*roots",
                              "14N" = "Control*' '*roots",
                              "13C + 15N" = "''^13*C*' '*'+'*''^15*N*' '*roots"), default = label_parsed)

# Setting up some functions
color_mapping <- function(P_level) {
  if (grepl("0", P_level)) {
    return(p0)  
  } else if (grepl("188", P_level)) {
    return(p188) 
  }
  
}

std <- function(x) sd(x)/sqrt(length(x))

# Creating a table with microbial biomass C data
metadata_exp1 <- data.frame(sample_data(exp1_ps_16s_r))

# MICROBIAL BIOMASS C
# Summarizing data
mbc_summ <- metadata_exp1 %>%
  filter(!is.na(MBC_mg.g)) %>%
  mutate(MBC_mg.g = MBC_mg.g/10) %>%
  group_by(Day, P.level) %>%
  summarise_at(vars(MBC_mg.g), funs(mean(.), std(.), sd(.), median(.), IQR(.)))
mbc_summ$color <- sapply(mbc_summ$P.level, color_mapping)


# Calculating difference between days
stats_mbc <- metadata_exp1 %>%
  filter(MBC_mg.g != "NA") %>%
  group_by(P.level) %>%
  wilcox_test(data = ., formula = MBC_mg.g ~ Day, ref.group = "2") %>%
  filter(p < 0.05)
stats_mbc$color <- sapply(stats_mbc$P.level, color_mapping)

# between treatments
metadata_exp1 %>%
  filter(MBC_mg.g != "NA") %>%
  group_by(Day) %>%
  wilcox_test(data = ., formula = MBC_mg.g ~ P.level, ref.group = "0") %>%
  filter(p < 0.05)

metadata_exp1 %>%
  wilcox_test(data = ., formula = MBC_mg.g ~ P.level, ref.group = "0", paired = FALSE) %>%
  add_significance()


# Building a plot
mbc_plot_exp1 <- ggplot(mbc_summ, aes(x = Day, y = mean, group = P.level, color = mbc_summ$color)) +
  geom_line(aes(color = mbc_summ$color, group = P.level)) + geom_point(aes(color = mbc_summ$color, group = P.level), size = 3) +
  geom_errorbar(aes(ymin = mean-std, ymax = mean+std,color = mbc_summ$color),
                width =.1, position = position_dodge(0)) +
  labs(x = "Time (d)", y = "MBC (µg/g)") + my_theme + theme(axis.ticks = element_blank()) +
  stat_pvalue_manual(stats_mbc, label = "p.adj.signif", hide.ns = TRUE, y.position = c(255, 290, 225), remove.bracket = T, size = 6, color = "color") +
  scale_color_identity(breaks = c(p188, p0), labels = c("188", "0"), guide = "legend") +
  guides(color = guide_legend(title = "SP (kg/ha)", override.aes = aes(label = ""))) +
  scale_x_discrete(labels = c("2", "3 *", "7")) +
  ylab(expression(bold("MBC (µg "*g^{-1}*")"))) +
  ylim(0, 320)

mbc_plot_exp1

# HWEC
# summarizing data for total organic carbon
tc_summ <- metadata_exp1 %>%
  filter(!is.na(TC)) %>%
  #mutate(TC = TC/1000) %>%
  group_by(Day, P.level) %>%
  summarise_at(vars(TC), funs(mean(.), std(.)))
tc_summ$color <- sapply(tc_summ$P.level, color_mapping)

metadata_exp1 %>% # no difference
  filter(TC != "NA") %>%
  group_by(Day) %>%
  wilcox_test(data = ., formula = TC ~ P.level, ref.group = "0") %>%
  filter(p < 0.05)

metadata_exp1 %>%
  wilcox_test(data = ., formula = TC ~ P.level, ref.group = "0", paired = FALSE) %>%
  add_significance()

# building plot for toc
tc_plot_exp1 <- ggplot(tc_summ, aes(x = Day, y = mean, group = P.level, color = P.level)) +
  geom_line(aes(color = tc_summ$color)) + geom_point(aes(color = tc_summ$color),size = 3) +
  geom_errorbar(data = tc_summ, aes(ymin = mean-std, ymax = mean+std,color = tc_summ$color),
                width =.1, position = position_dodge(0)) +
  labs(x = "Time (d)") + my_theme + theme(axis.ticks = element_blank()) +
  scale_color_identity(breaks = c(p188, p0), labels = c("188", "0"), guide = "legend") +
  guides(color = guide_legend(title = "SP (kg/ha)", override.aes = aes(label = ""))) +
  ylab(expression(bold("HWEC (µg "*g^{-1}*")")))
tc_plot_exp1

# GHG FLUX
flux_summ <- metadata_exp1 %>%
  filter(!is.na(C_CO2_hour)) %>%
  group_by(Day, P.level) %>%
  summarise_at(vars(C_CO2_hour, N_N2O_hour), funs(mean(.), std(.)))
flux_summ$color <- sapply(flux_summ$P.level, color_mapping)

# Calculating difference between days for CO2 & N2O
stats_co2 <- metadata_exp1 %>%
  group_by(P.level) %>%
  t_test(data = ., formula = C_CO2_hour ~ Day, ref.group = "1") %>%
  adjust_pvalue() %>%
  filter(p.adj < 0.05) %>%
  add_significance()
stats_co2$color <- sapply(stats_co2$P.level, color_mapping)

metadata_exp1 %>%
  t_test(data = ., formula = C_CO2_hour ~ P.level, ref.group = "0") %>%
  add_significance()


stats_n2o <- metadata_exp1 %>%
  group_by(P.level) %>%
  t_test(data = ., formula = N_N2O_hour ~ Day, ref.group = "1") %>%
  adjust_pvalue() %>%
  filter(p.adj < 0.05) %>%
  add_significance()
stats_n2o$color <- sapply(stats_n2o$P.level, color_mapping)


metadata_exp1 %>%
  t_test(data = ., formula = N_N2O_hour ~ P.level, ref.group = "0") %>%
  add_significance()

# building a plot for CO2

co2_plot_exp1 <- ggplot(flux_summ, aes(x = Day, y = C_CO2_hour_mean, group = P.level, color = P.level)) +
  geom_line(aes(color = flux_summ$color)) + geom_point(aes(color = flux_summ$color), size = 3) +
  geom_errorbar(aes(ymin = C_CO2_hour_mean - C_CO2_hour_std, ymax = C_CO2_hour_mean + C_CO2_hour_std, color = flux_summ$color), 
                width =.1, position = position_dodge(0.05)) +
  labs(x = "Time (d)", color = expression("SP (kg "*ha^{-1}*")")) +
  my_theme + theme(axis.ticks = element_blank()) +
  stat_pvalue_manual(stats_co2, label = "p.adj.signif", hide.ns = TRUE, y.position = c(1.6, 1.9, 1.2, 1.25, 1.35, 1.1), remove.bracket = T, size = 6, color = "color") +
  scale_color_identity(breaks = c(p188, p0), labels = c("188", "0"), guide = "legend") +
  guides(color = guide_legend(title = expression("SSP (kg "*ha^{-1}*' '*year^{-1}*")"), override.aes = aes(label = ""))) +
  ylim(0, 2.5) +
  ylab(expression(bold(CO[2]~"flux (µg C "*g^{-1}*h^{-1}*")")))
co2_plot_exp1


# ... for N2O
n2o_plot_exp1 <- ggplot(flux_summ, aes(x = Day, y = N_N2O_hour_mean, group = P.level, color = P.level)) +
  geom_line(aes(color = flux_summ$color)) + geom_point(aes(color = flux_summ$color), size = 2) +
  geom_errorbar(aes(ymin = N_N2O_hour_mean - N_N2O_hour_std, ymax = N_N2O_hour_mean + N_N2O_hour_std, color = flux_summ$color), 
                width =.1, position = position_dodge(0.05)) +
  labs(x = "Time (d)", color = "SP (kg/ha)") + 
  my_theme + theme(axis.ticks = element_blank()) +
  ylim(0, NA) +
  stat_pvalue_manual(stats_n2o, label = "p.adj.signif", hide.ns = TRUE, y.position = c(0.00026, 0.00037, 0.00035, 0.00012), remove.bracket = T, size = 6, color = "color") +
  scale_color_identity(breaks = c(p188, p0), labels = c("188", "0"), guide = "legend") +
  guides(color = guide_legend(title = "SP (kg/ha)", override.aes = aes(label = ""))) +
  ylab(expression(bold(N[2]*O~"flux (µg N "*g^{-1}*h^{-1}*")")))
n2o_plot_exp1

# INORGANIC N
# summarize 
kcl_summ <- metadata_exp1 %>%
  filter(!is.na(NH4_N)) %>%
  group_by(Day, P.level) %>%
  summarise_at(vars(NH4_N, NO3_N), funs(mean(.), std(.)))
kcl_summ$color <- sapply(kcl_summ$P.level, color_mapping)

# calculating difference between days for NH4 N
stats_nh4 <- metadata_exp1 %>%
  group_by(P.level) %>%
  t_test(data = ., formula = NH4_N ~ Day, ref.group = "1", p.adjust.method = "none") %>%
  filter(p < 0.05) %>%
  add_significance()
stats_nh4$color <- sapply(stats_nh4$P.level, color_mapping)

# ... and NO3 N
stats_no3 <- metadata_exp1 %>%
  group_by(P.level) %>%
  t_test(data = ., formula = NO3_N ~ Day, ref.group = "1", p.adjust.method = "none") %>%
  filter(p < 0.05) %>%
  add_significance()
stats_no3$color <- sapply(stats_no3$P.level, color_mapping)


metadata_exp1 %>%
  group_by(Day) %>%
  t_test(data = ., formula = NH4_N ~ P.level, ref.group = "0", p.adjust.method = "none") %>%
  add_significance()

# plotting NH4 N
nh4_plot_exp1 <- ggplot(kcl_summ, aes(x = Day, y = NH4_N_mean, group = P.level, color = P.level)) +
  geom_line() + geom_point(size = 3) +
  geom_errorbar(aes(ymin = NH4_N_mean - NH4_N_std, ymax = NH4_N_mean + NH4_N_std), 
                width =.1, position = position_dodge(0.05)) +
  labs(x = "Time (d)") +
  my_theme + theme(axis.ticks = element_blank()) + 
  ylim(20, 35) +
  scale_color_manual(name = "SP (kg/ha)", values = c(p0, p188)) +
  ggtitle("B") +
  ylab(expression(bold(NH[4]*' - '*N*' '*"(µg "*g^{-1}*")"))) +
  scale_x_discrete(labels = c("1", "2 *", "3", "4 **", "7")) +
  stat_pvalue_manual(stats_nh4, label = "p.adj.signif", hide.ns = TRUE, y.position = 24, remove.bracket = T, size = 6, color = "#708A58") 
nh4_plot_exp1

# ... and NO3 N

metadata_exp1 %>%
  group_by(Day) %>%
  t_test(data = ., formula = NO3_N ~ P.level, ref.group = "0") %>%
  add_significance()

no3_plot_exp1 <- ggplot(kcl_summ, aes(x = Day, y = NO3_N_mean, group = P.level, color = P.level)) +
  geom_line() + geom_point(size = 3) +
  geom_errorbar(aes(ymin = NO3_N_mean - NO3_N_std, ymax = NO3_N_mean + NO3_N_std), 
                width =.1, position = position_dodge(0.05)) +
  labs(x = "Time (d)") +
  my_theme + theme(axis.ticks = element_blank()) + 
  stat_pvalue_manual(stats_no3, label = "p.adj.signif", hide.ns = TRUE, y.position = 70, remove.bracket = T, size = 6, color = "#708A58") +
  scale_color_manual(name = "SP (kg/ha)", values = c(p0, p188)) +
  ylab(expression(bold(NO[3]*' - '*N*' '*"(µg "*g^{-1}*")")))
no3_plot_exp1

# 15 N N2O
n15_summ <- metadata_exp1 %>%
  filter(!is.na(X15N_atom.)) %>%
  group_by(Day, P.level, N.treatment, OM.trt) %>%
  summarise_at(vars(X15N_atom.), funs(median(.), IQR(.), mean(.), std(.)))
n15_summ$color = sapply(n15_summ$P.level, color_mapping)

# Calculating difference between days for 15 N
stats_n15 <- metadata_exp1 %>%
  group_by(P.level, N.treatment, OM.trt) %>%
  t_test(data = ., formula = X15N_atom. ~ Day, ref.group = "1", p.adjust.method = "none") %>%
  filter(p < 0.05) %>%
  add_significance()
stats_n15$color <- sapply(stats_n15$P.level, color_mapping)

metadata_exp1 %>% # no difference
  filter(X15N_atom. != "NA") %>%
  group_by(Day) %>%
  wilcox_test(data = ., formula = X15N_atom. ~ P.level, ref.group = "0", p.adjust.method = "none") %>%
  filter(p < 0.05)

metadata_exp1 %>%
  wilcox_test(data = ., formula = X15N_atom. ~ P.level, ref.group = "0", paired = FALSE) %>%
  add_significance()

# Building a plot for 15 N
library(ggh4x)
n15_plot_exp1 <- ggplot(n15_summ, aes(x = Day, y = mean, group = P.level, color = color)) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - std, ymax = mean + std),
                width = 0.1, position = position_dodge(0.05)) +
  my_theme +
  facet_wrap(~ N.treatment, labeller = labeller(N.treatment = Ntreatments), scales = "free_y") +
  stat_pvalue_manual(data = stats_n15, label = "p.adj.signif", hide.ns = TRUE,
                     y.position = c(0.52, 1.1, 0.9, 0.6, 1),
                     remove.bracket = TRUE, size = 6, color = "color") +
  scale_color_identity(breaks = c(p188, p0),
                       labels = c("188", "0"),
                       guide = "legend") +
  guides(color = guide_legend(title = "SP (kg/ha)", override.aes = aes(label = ""))) +
  ylab(expression(bold(N[2]*O~' '^{"15"}~N*' '*"(atom%)"))) +
  xlab("Time (d)") +
  theme(axis.ticks = element_blank()) +
  ggh4x::facetted_pos_scales(
    y = list(
      N.treatment == "14NH414NO3" ~ scale_y_continuous(limits = c(0, 1.5)),
      N.treatment == "14NH415NO3" ~ scale_y_continuous(limits = c(0, 8)),
      N.treatment == "15NH414NO3" ~ scale_y_continuous(limits = c(0, 1.5))
    )
  )
