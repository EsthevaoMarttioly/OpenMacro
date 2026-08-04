#=
# ---------------------------------------------------------------------------
# Results: prints every table, writes them to output/, then draws the figures.
# ---------------------------------------------------------------------------
#=

library(ggplot2)


# ---------------------------------------------------------------------------
cat("\n=== Q2: Brazil,", min(br$year), "-", max(br$year), "===\n")
print(as.data.frame(q2), digits = 3)
cat("\nMaddison: both countries annual from", start, "to", max(yrs), "\n")
cat("sigma_BR/sigma_US:", min(ratio$year), "-", max(ratio$year), "| mean/min/max",
    round(c(mean(ratio$ratio), min(ratio$ratio), max(ratio$ratio)), 2), "\n")

cat("\n=== Q3: sigma_C / sigma_Y ===\n")
print(as.data.frame(q3), digits = 3)
print(as.data.frame(q3_windows), digits = 3)

cat("\n=== Q4:", nrow(country), "countries,", er_from, "-", er_to, "===\n")
print(as.data.frame(q4), digits = 3)
print(as.data.frame(q4_income), digits = 3)

cat("\n=== Q5: dy_t = rho * dy_{t-1} + e_t ===\n")
print(as.data.frame(q5), digits = 3)
cat("implied_ratio is (1+r)/(1+r-rho) with r =", r_ss,
    "- compare with sigma_C/sigma_Y above.\n\n")

write.csv(q2,        "output/q2_tb_correlations.csv",       row.names = FALSE)
write.csv(ratio,     "output/q2_sigma_br_us.csv",           row.names = FALSE)
write.csv(rolling,   "output/q3_sigma_ratio_rolling.csv",   row.names = FALSE)
write.csv(country,   "output/q4_er_openness.csv",           row.names = FALSE)
write.csv(q4_income, "output/q4_er_openness_by_income.csv", row.names = FALSE)
write.csv(q5,        "output/q5_ar1.csv",                   row.names = FALSE)

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

save_fig = function(g, name)
  ggsave(paste0("output/", name, ".png"), g, width = 8, height = 5, dpi = 150, bg = "white")


# ---- Graphics -------------------------------------------------------------

# Q2
g = br %>% filter(!is.na(cy)) %>%
  ggplot(aes(year)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(aes(y = 100 * cy, colour = "Output gap (Hamilton)"), linewidth = 1) +
  geom_line(aes(y = tby, colour = "tb/y (raw)"), linewidth = 1) +
  scale_colour_manual(values = c(`Output gap (Hamilton)` = "#1f4e79",
                                 `tb/y (raw)` = "#c0504d")) + mytheme +
  labs(title = "Brazil: output gap and the trade balance",
       subtitle = paste0("Only y is detrended. rho = ",
                         round(q2$`rho(tb/y, y)`, 2)),
       x = NULL, y = "% / p.p.", colour = "")
save_fig(g, "q2_gap_vs_tb")

g = ggplot(ratio, aes(year, ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(linewidth = 1.1, colour = "#2c8c62") + mytheme +
  labs(title = "Gap volatility: Brazil relative to the United States",
       subtitle = paste0("Rolling ", win, "-year sd of the Hamilton gap in log real ",
                         "GDP per capita, Maddison Project"),
       x = NULL, y = expression(sigma[y]^{BR} / sigma[y]^{US}))
save_fig(g, "q2_sigma_br_us")

g = gaps %>% filter(!is.na(sigma)) %>%
  ggplot(aes(year, 100 * sigma, colour = iso)) +
  geom_line(linewidth = 1.1) + mytheme +
  scale_colour_manual(values = c(BRA = "#2c8c62", USA = "#1f4e79")) +
  labs(title = "Output gap volatility over the long run",
       subtitle = paste0("Rolling ", win, "-year sd of the Hamilton gap, Maddison Project"),
       x = NULL, y = "%", colour = "")
save_fig(g, "q2_long_run_volatility")


# Q3
g = ggplot(rolling, aes(year, ratio, colour = window)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(linewidth = 1.1) + mytheme +
  scale_colour_manual(values = c("20 years" = "#2c8c62", "30 years" = "#1f4e79")) +
  labs(title = "Brazil: rolling sigma_C / sigma_Y",
       subtitle = paste0("Real per capita, Hamilton (2018) h = 2, p = 2, ",
                         "right-aligned window. Last observation: ", max(br$year), "."),
       x = NULL, y = expression(sigma[C] / sigma[Y]), colour = "Window")
save_fig(g, "q3_sigma_ratio_rolling")

g = br %>% filter(!is.na(cy)) %>%
  ggplot(aes(year)) + mytheme +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(aes(y = 100 * cy, colour = "Output"), linewidth = 1) +
  geom_line(aes(y = 100 * cc, colour = "Consumption"), linewidth = 1) +
  scale_colour_manual(values = c(Output = "#1f4e79", Consumption = "#c0504d")) +
  labs(title = "Brazil: cyclical components, real per capita",
       subtitle = "Hamilton (2018), h = 2, p = 2",
       x = NULL, y = "% deviation from trend", colour = "")
save_fig(g, "q3_brazil_cycles")


# Q4
g = ggplot(country, aes(openness, 100 * sigma_er)) +
  geom_point(aes(colour = income), alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.9) +
  geom_text(data = filter(country, iso %in% c("BRA", "ARG", "JPN", "CHE", "PAN",
                                              "SGP", "ZAF", "MEX", "TUR", "CHN",
                                              "IND", "DEU", "NLD", "HKG")),
            aes(label = iso), vjust = -0.9, size = 2.8, colour = "grey20") +
  scale_x_log10() + scale_y_log10() + mytheme +
  labs(title = "Exchange rate volatility and openness",
       subtitle = paste0("One point per country, ", er_from, "-", er_to,
                         ". Both axes on a log scale."),
       x = "Trade openness, (X + M) / Y (%)", colour = "",
       y = expression(sigma[ER] * ", sd of annual " * Delta * "log(LCU/US$) (%)"))
save_fig(g, "q4_er_vs_openness")

g = country %>%
  mutate(bucket = cut(openness, breaks = quantile(openness, 0:4 / 4),
                      labels = c("Most closed", "2nd", "3rd", "Most open"), include.lowest = TRUE)) %>%
  ggplot(aes(bucket, 100 * sigma_er)) + scale_y_log10() + mytheme +
  geom_boxplot(fill = "#1f4e79", alpha = 0.25, outlier.size = 1) +
  labs(title = "Exchange rate volatility by openness quartile",
       subtitle = paste0("Countries sorted by mean (X + M) / Y, ", er_from, "-", er_to),
       x = NULL, y = expression(sigma[ER] * " (%, log scale)"))
save_fig(g, "q4_er_by_openness_quartile")


# Q5
g = ggplot(acf_df, aes(lag, acf, colour = sample)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(linewidth = 1.1) + geom_point(size = 2) +
  scale_x_continuous(breaks = 1:8) + mytheme +
  scale_colour_manual(values = c(Annual = "#1f4e79", `Quarterly SA` = "#c0504d")) +
  labs(title = "Brazil: autocorrelation of GDP growth",
       subtitle = "An AR(1) implies a geometric decay from the first bar",
       x = "Order", y = "Autocorrelation", colour = "")
save_fig(g, "q5_growth_acf")

g = qtr %>%
  mutate(date = as.Date(paste0(year, "-", 3 * quarter - 2, "-01"))) %>%
  ggplot(aes(date, index)) +
  geom_line(linewidth = 1, colour = "#1f4e79") + mytheme +
  labs(title = "Brazil: quarterly GDP volume index", x = NULL, y = "Index",
       subtitle = "IBGE via BCB SGS 22109, seasonally adjusted, 1995 average = 100")
save_fig(g, "q5_brazil_gdp_quarterly")
