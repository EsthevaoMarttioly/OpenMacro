#=
# ---------------------------------------------------------------------------
# Results: writes the tables into output/tables/ and the figures into
# output/figures/. Nothing is printed to the console.
# ---------------------------------------------------------------------------
#=

library(ggplot2)
library(xtable)

# Whole float, so tasks.tex only needs \input{tables/<name>}.
# [H] from the float package pins the table where it is placed. fit = "width"
# shrinks a wide table, fit = "page" gives it a page of its own.
# Column names go through untouched, so they can carry math.
save_tex = function(d, name, caption, label, digits = 3,
                    fit = c("none", "width", "page")) {
  fit = match.arg(fit)
  tab = capture.output(
    print(xtable(as.data.frame(d), digits = digits),
          include.rownames = FALSE, floating = FALSE, booktabs = TRUE,
          comment = FALSE, sanitize.colnames.function = identity))

  body = switch(fit,
    none  = tab,
    width = c("\\resizebox{\\textwidth}{!}{%", tab, "}"),
    page  = c("\\resizebox*{!}{\\dimexpr\\textheight-2\\baselineskip\\relax}{%",
              tab, "}"))

  writeLines(c(paste0("\\begin{table}", if (fit == "page") "[p]" else "[H]"),
               "\\centering", body,
               paste0("\\caption{", caption, "}"),
               paste0("\\label{tab:", label, "}"),
               "\\end{table}"),
             paste0("output/tables/", name, ".tex"))
}


etable(q1, tex = TRUE, file = "output/tables/q1_debt_interest.tex", replace = TRUE,
       digits = 3, fitstat = ~ n + r2, placement = "H",
       headers = list(Panel = c("All", "All", "Developed", "Emerging")),
       title = "Regression with fixed effects of external debt on interest rate",
       label = "tab:q1",
       dict = c(interest = "Interest Rate (\\%)",
                `lag(debt)` = "$debt_{t-1}$ (\\% GDP)", iso = "Country"))

save_tex(setNames(q3_windows, c("Window", "From", "To", "Min", "Max", "Range")),
         "q3_windows",
         "Rolling $\\sigma_c/\\sigma_y$, Brazil.", "q3")

save_tex(setNames(q4_income, c("Income Group", "N", "Spearman",
                               "Median $\\sigma_{ER}$", "Median Openness")),
         "q4_income",
         "Correlation of $\\sigma_{ER}$ with openness, within income group.", "q4")

save_tex(setNames(q5, c("Sample", "$\\rho$", "P-value", "$\\sigma_\\varepsilon$",
                        "N", "Implied Ratio")),
         "q5_ar1",
         paste0("$\\Delta y_t = \\alpha + \\rho \\Delta y_{t-1} + ",
                "\\varepsilon_t$. Implied ratio uses $r = 0.04$."),
         "q5", digits = 4)

q6_names = c("Group / Means", "Trade Balance", "Interest Rate", "Debt",
             "Debt Implied", "Ratio", "From", "To")
q6_caption = paste("Ratio of trade balance and interest rate times external",
                   "debt, in average. Trade balance and debt in \\% of GDP.")

save_tex(setNames(q6_regions, q6_names), "q6_weo_regions",
         q6_caption, "ratio", fit = "width")

save_tex(setNames(q6_countries, sub("Group", "Country", q6_names)), "q6_countries",
         paste(q6_caption, "By country."), "ratio_wb", fit = "page")

# Numbers quoted in the prose, as macros, so the text follows the data
writeLines(c(
  sprintf("\\newcommand{\\qtwoRho}{%.3f}",       q2$`rho(tb/y, y)`),
  sprintf("\\newcommand{\\qtwoFrom}{%d}",        min(br$year)),
  sprintf("\\newcommand{\\qtwoTo}{%d}",          max(br$year)),
  sprintf("\\newcommand{\\qtwoRatioMean}{%.2f}", mean(ratio$ratio)),
  sprintf("\\newcommand{\\qtwoRatioMin}{%.2f}",  min(ratio$ratio)),
  sprintf("\\newcommand{\\qtwoRatioMax}{%.2f}",  max(ratio$ratio)),
  sprintf("\\newcommand{\\qfourSlope}{%.3f}",    q4$slope),
  sprintf("\\newcommand{\\qfourRtwo}{%.3f}",     q4$r2),
  sprintf("\\newcommand{\\qfourSpearman}{%.3f}", q4$spearman),
  sprintf("\\newcommand{\\qfourN}{%d}",          q4$n)
), "output/tables/numbers.tex")

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

save_fig = function(g, name, height = 5)
  ggsave(paste0("output/figures/", name, ".png"), g, width = 8, height = height,
         dpi = 150, bg = "white")


# ---- Graphics -------------------------------------------------------------

# Q1
g = q1_country %>%
  mutate(name = fct_reorder(name, beta1)) %>%
  ggplot(aes(beta1, name)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = beta1 - 1.96 * se, xmax = beta1 + 1.96 * se),
                  linewidth = 0.4, size = 0.25, colour = "#1f4e79") +
  mytheme +
  labs(title = "Interest rate on lagged external debt, country by country",
       subtitle = "OLS per country, 95% confidence interval",
       x = expression(hat(beta)[1] * ", p.p. of rate per p.p. of GDP"), y = NULL)
save_fig(g, "q1_by_country", height = 9)


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


# Q6
g = q6_countries %>%
  filter(ratio > 0, is.finite(ratio)) %>%
  ggplot(aes(ratio)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_histogram(bins = 40, fill = "#1f4e79", alpha = 0.7) +
  scale_x_log10() + mytheme +
  labs(title = "Steady-state debt implied by tb = r d, against the actual mean",
       subtitle = "Both in % of GDP. Ratio of implied to actual, log scale. Negative ratios dropped.",
       x = expression(bar(d)[implied] / bar(d)[actual]), y = "Countries")
save_fig(g, "q6_ratio_hist")
