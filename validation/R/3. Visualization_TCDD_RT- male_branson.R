########
## Plot simulations
###

# ridgeplot for all tisssues

library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)

# 1. Filtrer kun tid og konsentrasjoner
sim_data_long <- sim_data_branson %>%
  select(time, starts_with("C_")) %>%
  pivot_longer(-time, names_to = "Tissue", values_to = "Concentration")

# 2. Fjern total, vann og urin hvis ønskelig
sim_data_long <- sim_data_long %>%
  filter(!Tissue %in% c("C_tot", "C_water", "C_urine", "C_plasma"))

# 3. Formatér tissue-navn for lesbarhet (valgfritt)
sim_data_long$Tissue <- gsub("C_", "", sim_data_long$Tissue)
sim_data_long$Tissue <- factor(sim_data_long$Tissue, levels = unique(sim_data_long$Tissue))

# 4. Lag ridgeplot
ridge_plot <- ggplot(sim_data_long, aes(x = time, y = Tissue, height = Concentration, group = Tissue, fill = Tissue)) +
  geom_ridgeline(scale = 1, alpha = 0.8) +
  labs(
    title = "Time profile of chemical concentration by tissue",
    x = "Time",
    y = "Tissue"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


print(ridge_plot)

#####
# Facet wrap plot

# 1. Data prep
sim_data_long <- sim_data_branson %>%
  select(time, starts_with("C_")) %>%
  pivot_longer(cols = -time, names_to = "Tissue", values_to = "Concentration") %>%
  filter(!Tissue %in% c("C_tot", "C_water", "C_urine", "C_plasma")) %>%
  mutate(
    Tissue = gsub("C_", "", Tissue),
    Tissue = factor(Tissue)
  )

# 2. Clean and aesthetic lineplot layout
line_plot <- ggplot(sim_data_long, aes(x = time, y = Concentration)) +
  geom_line(color = "#1f77b4", size = 0.8) +
  facet_wrap(~ Tissue, scales = "fixed", ncol = 2) + # max scale use: scales = "free_y" , same Y: "fixed"
  labs(
    title = "Tissue-specific chemical concentration over time",
    x = "Time (days)",
    y = "Concentration (\u00b5g/g)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 9)
  )

print(line_plot)

############
# Colored overlay - linear

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsci)

# 1. Data prep
sim_data_long <- sim_data_branson %>%
  select(time, starts_with("C_")) %>%
  pivot_longer(cols = -time, names_to = "Tissue", values_to = "Concentration") %>%
  filter(!Tissue %in% c("C_tot", "C_water", "C_urine", "C_plasma")) %>%
  mutate(
    Tissue = gsub("C_", "", Tissue),
    Tissue = factor(Tissue)
  )

# 2. Overlagt linjediagram med farger per vev
combined_plot <- ggplot(sim_data_long, aes(x = time, y = Concentration, color = Tissue)) +
  geom_line(size = 1) +
  labs(
    title = "Chemical concentration in tissues over time",
    x = "Time (days)",
    y = "Concentration (\u00b5g/g)",
    color = "Tissue"
  ) +
  # scale_x_log10() +
  scale_color_d3("category20") +  # max 20 colors
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

print(combined_plot)
