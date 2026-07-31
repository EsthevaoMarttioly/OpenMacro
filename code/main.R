#=
# ---------------------------------------------------------------------------
# DESCRIPTION
# 
# ---------------------------------------------------------------------------
#=

# ---- Packages -------------------------------------------------------------
# renv::restore()     # run once to install the locked versions ('1' or 'Y')
library(tidyverse)
library(readxl)


# ---- Config ---------------------------------------------------------------



# ---------------------------------------------------------------------------
# 1. 





# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
mytheme = theme(legend.position = "bottom",
                plot.title = element_text(size = 12, face = "bold"),
                plot.subtitle = element_text(size = 10),
                panel.background = element_rect(fill = "transparent", colour = "black",
                                                linewidth = 0.5, linetype = "solid"),
                panel.grid.major.y = element_line(colour = "grey", linewidth = 0.5),
                panel.grid.minor.y = element_line(colour = "grey", linewidth = 0.5),
                panel.grid = element_line(colour = "grey98"),
                panel.grid.major.x = element_line(colour = "transparent"),
                panel.grid.minor.x = element_line(colour = "transparent"),
                axis.text = element_text(colour = "black", size = 9),
                strip.background = element_rect(fill = "grey95", colour = "black"),
                strip.text = element_text(colour = "black", size = 9))

save_fig = function(g, name) ggsave(paste0("output/figures/", name, ".png"),
                                    g, width = 8, height = 5, dpi = 150)


# ---- Graphics -------------------------------------------------------------
g = ggplot(wage_dist, aes(wage, density, colour = sector)) +
  geom_line(linewidth = 1.2) + mytheme +
  geom_line(aes(y = lognormal), linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = 1518, color = "black", linetype = "dashed", linewidth = 0.8) +
  scale_colour_manual(values = col_sector) +
  coord_cartesian(xlim = c(0, 7e3)) +
  labs(title = paste("PNAD", max(year), "- Wage Distribution vs Log-Normal"),
       x = "Monthly Wage (R$)", y = "Density", colour = "")
save_fig(g, "wage_distribution")



