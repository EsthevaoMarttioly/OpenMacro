#=
# ---------------------------------------------------------------------------
# DESCRIPTION
# Open Macro mini-tasks. Q7 to Q9 are in dynare/.
# The data lives in code/data.R and only has to run when you want to refresh it. 
# ---------------------------------------------------------------------------
#=

# ---- Packages -------------------------------------------------------------
# renv::restore()     # run once to install the locked versions ('1' or 'Y')
library(tidyverse)
library(fixest)


# ---- Config ---------------------------------------------------------------
if (!file.exists("data/brazil_annual.csv")) source("code/data.R")

win     = 30
wins    = c(20, 30)
er_from = 1990         # Q4, window for the cross-country moments
er_to   = 2024
r_ss    = 0.04         # Q5, steady-state real rate for the implied sigma ratio
q5_from = 1996         # Q5, annual sample starts with the quarterly one

hamilton_cycle = function(x, h = 2, p = 2) {
  n = length(x)
  y = x[(h + p):n]
  X = sapply(0:(p - 1), function(j) x[(h + p - 1 - j):(n - 1 - j)])
  c(rep(NA_real_, h + p - 1), residuals(lm(y ~ X)))
}

roll_sd = function(x, win) {      # Rolling window for Q3
  out = rep(NA_real_, length(x))
  for (t in win:length(x)) out[t] = sd(x[(t - win + 1):t])
  out
}

gmean_pct = function(x) 100 * (exp(mean(log1p(x / 100), na.rm = TRUE)) - 1)

# Databases
br = read.csv("data/brazil_annual.csv") %>%
  filter(!is.na(gdp_pc), !is.na(cons_pc), !is.na(tby)) %>%
  mutate(cy = hamilton_cycle(log(gdp_pc)),     # output gap
         cc = hamilton_cycle(log(cons_pc)))    # consumption gap, needed for Q3

mpd = read.csv("data/maddison_gdp_pc.csv")

panel = read.csv("data/country_panel.csv") %>%
  filter(between(year, er_from, er_to), iso != "USA")   # US is the numeraire

qtr = read.csv("data/brazil_gdp_quarterly_sa.csv")

debt = read.csv("data/debt_panel.csv") %>%
  filter(!iso %in% c("UKR", "SRB", "BIH", "CHN", "GEO"))  # countries with broken debt series

weo = read.csv("data/weo_regions.csv")


# ---------------------------------------------------------------------------
# Q1. Interest Rates and External Debt
q1_reg = debt %>%
  filter(!is.na(interest), !is.na(debt)) %>%
  mutate(group = if_else(income == "High income", "Developed", "Emerging"))

q1 = list(
  `Pooled OLS, all`        = feols(interest ~ lag(debt), q1_reg, vcov = ~ iso),
  `Country FE, all`        = feols(interest ~ lag(debt) | iso, q1_reg, vcov = ~ iso),
  `Country FE, developed`  = feols(interest ~ lag(debt) | iso,
                                   filter(q1_reg, group == "Developed"), vcov = ~ iso),
  `Country FE, emerging`   = feols(interest ~ lag(debt) | iso,
                                   filter(q1_reg, group == "Emerging"), vcov = ~ iso)
)

# One regression per country, to see how much the pooled slope hides
q1_country = q1_reg %>% group_by(iso, name) %>%
  filter(n() >= 10) %>%
  reframe(broom::tidy(lm(interest ~ lag(debt))) %>% slice(2),
          from = min(year), to = max(year)) %>%
  select(iso, name, beta1 = estimate, se = std.error, t = statistic, from, to) %>%
  arrange(beta1)


# ---------------------------------------------------------------------------
# Q2. GDP and Trade Balance
q2 = tibble(`rho(tb/y, y)` = cor(br$tby, br$cy, use = "complete.obs"))

# Extra: Long series for Brazilian GDP
yrs = mpd %>% count(year) %>% filter(n == 2) %>% pull(year) %>% sort()
start = yrs[max(which(c(2, diff(yrs)) > 1))]

gaps = mpd %>% filter(year %in% yrs, year >= start) %>%
  group_by(iso) %>% arrange(year, .by_group = TRUE) %>%
  mutate(cy    = hamilton_cycle(log(gdp_pc)),
         sigma = roll_sd(cy, win)) %>% ungroup()

ratio = gaps %>% select(iso, year, sigma) %>%
  pivot_wider(names_from = iso, values_from = sigma) %>%
  filter(!is.na(BRA), !is.na(USA)) %>% mutate(ratio = BRA / USA)


# ---------------------------------------------------------------------------
# Q3. Sigma Ratios
q3 = tibble(`sigma_Y (%)`     = 100 * sd(br$cy, na.rm = TRUE),
            `sigma_C (%)`     = 100 * sd(br$cc, na.rm = TRUE),
            `sigma_C/sigma_Y` = sd(br$cc, na.rm = TRUE) / sd(br$cy, na.rm = TRUE),
            `rho(c, y)`       = cor(br$cc, br$cy, use = "complete.obs"))

rolling = map_dfr(wins, function(w)
  transmute(br, year, window = paste0(w, " years"),
                      ratio  = roll_sd(cc, w) / roll_sd(cy, w))) %>% filter(!is.na(ratio))

# how much of the 20-year swing survives the longer window
q3_windows = rolling %>%
  group_by(window) %>%
  summarise(from = min(year), to = max(year),
            min = min(ratio), max = max(ratio), range = max(ratio) - min(ratio))


# ---------------------------------------------------------------------------
# Q4. Exchange Rate Volatility vs Trade Volume (optional)
#     sigma_ER is heavily right-skewed, so the regression is log-log.
country = panel %>% group_by(iso, name, region, income) %>%
  summarise(sigma_er = sd(if_else(year - lag(year) == 1,
                                  log(er) - lag(log(er)), NA_real_), na.rm = TRUE),
            openness = mean(open, na.rm = TRUE),
            n = sum(!is.na(er) & !is.na(open)), .groups = "drop") %>%
  filter(n >= 20, is.finite(sigma_er), is.finite(openness), sigma_er > 0, openness > 0)

q4_fit = lm(log(sigma_er) ~ log(openness), data = country)

q4 = tibble(slope    = coef(q4_fit)[2],
            se       = summary(q4_fit)$coefficients[2, 2],
            r2       = summary(q4_fit)$r.squared,
            pearson  = cor(log(country$sigma_er), log(country$openness)),
            spearman = cor(country$sigma_er, country$openness, method = "spearman"),
            n        = nrow(country))

# same thing within income group, to check it is not just rich vs poor
q4_income = country %>% group_by(income) %>% filter(n() >= 10) %>%
  summarise(n = n(),
            spearman     = cor(sigma_er, openness, method = "spearman"),
            median_sigma = 100 * median(sigma_er),
            median_open  = median(openness), .groups = "drop")


# ---------------------------------------------------------------------------
# Q5. GDP as an AR(1)
growth = list(Annual = br$dy[br$year >= q5_from], `Quarterly SA` = qtr$dy)

q5 = imap_dfr(growth, function(dy, label) {
  m = lm(dy ~ dplyr::lag(dy))
  broom::tidy(m) %>% slice(2) %>%
    transmute(sample = label, rho = estimate, p = p.value,
              sigma_e = sigma(m), n = nobs(m),
              implied_ratio = (1 + r_ss) / (1 + r_ss - rho))})

acf_df = imap_dfr(growth, ~ tibble(sample = .y, lag = 1:8,
  acf = acf(.x, lag.max = 8, plot = FALSE, na.action = na.pass)$acf[-1]))


# ---------------------------------------------------------------------------
# Q6. External Debt Calibration
# tb and debt are both in % of GDP, so d_implied is too
q6_summary = function(d) {
  d %>% summarise(tb_bar   = mean(tb, na.rm = TRUE),
                  r_bar    = gmean_pct(interest),
                  debt_bar = mean(debt, na.rm = TRUE),
                  from     = min(year), to = max(year), .groups = "drop") %>%
    mutate(d_implied = tb_bar / (r_bar / 100),   # tb_bar = r * d_bar
           ratio     = d_implied / debt_bar) %>%
    select(unit_name, tb_bar, r_bar, debt_bar, d_implied, ratio, from, to)
}

q6_regions = weo %>% group_by(unit_name = region) %>% q6_summary()

q6_countries = debt %>% filter(!is.na(interest)) %>%
  group_by(unit_name = name) %>% q6_summary() %>% arrange(ratio)


# ---- Results --------------------------------------------------------------
source("code/results.R")

