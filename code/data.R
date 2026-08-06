#=
# ---------------------------------------------------------------------------
#  brazil_gdp_quarterly.csv     volume index, NSA                 1996+   SIDRA
#  brazil_gdp_quarterly_sa.csv  volume index, seasonally adj.     1996+   BCB
#  brazil_annual.csv            GDP, C, pop, tb/y                 1960+   WDI
#  country_panel.csv            exchange rate, openness           1960+   WDI
#  maddison_gdp_pc.csv          real GDP per capita, BRA and USA  1800+   MPD
#  debt_panel.csv               tb, external debt, interest       1970+   WB/OECD/IMF
#  weo_regions.csv              tb, external debt, interest       1980+   IMF WEO
#
# The IMF extracts and Debt_Interest.csv are hand-downloaded. See the README.
# ---------------------------------------------------------------------------
#=

library(tidyverse)
library(sidrar)
library(rbcb)
library(WDI)


# ---------------------------------------------------------------------------
# SIDRA / IBGE: quarterly GDP
sidra_raw = get_sidra(api = "/t/1620/n1/all/v/583/p/all/c11255/90707")
period_col = grep("Trimestre \\(C", names(sidra_raw), value = TRUE)[1]

gdp_q = tibble(period = as.character(sidra_raw[[period_col]]),
               index  = as.numeric(sidra_raw$Valor)) %>%
  mutate(year    = as.integer(substr(period, 1, 4)),
         quarter = as.integer(substr(period, 5, 6))) %>%
  filter(!is.na(index)) %>%
  select(period, year, quarter, index) %>% arrange(year, quarter)

write.csv(gdp_q, "data/brazil_gdp_quarterly.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# BCB: IBGE quarterly GDP, but seasonally adjusted.
gdp_q_sa = get_series(c(index = 22109)) %>% as_tibble() %>%
  transmute(year    = as.integer(format(date, "%Y")),
            quarter = (as.integer(format(date, "%m")) + 2) / 3, index) %>%
  arrange(year, quarter) %>% mutate(dy = c(NA, diff(log(index))))

write.csv(gdp_q_sa, "data/brazil_gdp_quarterly_sa.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# World Bank
# https://api.worldbank.org/v2/country/BRA/indicator/NY.GDP.MKTP.KN?format=json
options(timeout = 900)

wb_ind = c(gdp  = "NY.GDP.MKTP.KN",  # GDP, constant LCU
           cons = "NE.CON.PRVT.KN",  # household consumption, constant LCU
           pop  = "SP.POP.TOTL",
           tby  = "NE.RSB.GNFS.ZS",  # external balance G&S, % of GDP
           er   = "PA.NUS.FCRF",     # LCU per US$, period average
           open = "NE.TRD.GNFS.ZS",  # (X + M) / Y, %
           tb   = "BN.GSR.GNFS.CD",  # net trade in G&S, current US$
           debt = "DT.DOD.DECT.CD",  # external debt stocks, current US$
           gdp_usd = "NY.GDP.MKTP.CD")  # GDP, current US$, to scale the two above

wb_get = function(name, code, tries = 5) {
  for (i in seq_len(tries)) {
    x = try(WDI(country = "all", start = 1960, indicator = setNames(code, name)),
            silent = TRUE)
    if (!inherits(x, "try-error") && name %in% names(x))
      return(as_tibble(x) %>% select(iso = iso3c, year, all_of(name)))
    message("  ", code, " dropped by the server, retry ", i, "/", tries)
    Sys.sleep(5)
  }
  stop("World Bank never served ", code, ". Try again in a few minutes.")
}

wb = imap(wb_ind, ~ wb_get(.y, .x)) %>% reduce(full_join, by = c("iso", "year"))

# Region and income come from the table cached inside the WDI package
meta = as_tibble(WDI_data$country) %>%
  filter(region != "Aggregates") %>%
  select(iso = iso3c, name = country, region, income)

brazil = wb %>% filter(iso == "BRA") %>%
  select(year, gdp, cons, pop, tby) %>%
  arrange(year) %>% mutate(gdp_pc  = gdp / pop,
                           cons_pc = cons / pop,
                           dy      = c(NA, diff(log(gdp))))  # growth, for Q5

# Cross-country: Exchange Rates and Trade Volume
panel = wb %>% inner_join(meta, by = "iso") %>%
  select(iso, name, region, income, year, er, open) %>% arrange(iso, year)

write.csv(brazil, "data/brazil_annual.csv", row.names = FALSE)
write.csv(panel,  "data/country_panel.csv", row.names = FALSE)


# Average interest on new external debt commitments
ids_file = list.files("data", "Debt_Interest", full.names = TRUE)[1]
if (is.na(ids_file)) stop("Debt_Interest.csv missing from data/. See the README.")

r_wb = read.csv(ids_file, check.names = FALSE, colClasses = "character") %>%
  filter(`Series Code` == "DT.INR.DPPG") %>%
  select(iso = `Country Code`, matches("^[0-9]{4} ")) %>%
  pivot_longer(-iso, names_to = "year", values_to = "r_wb") %>%
  mutate(year = as.integer(substr(year, 1, 4)),
         r_wb = suppressWarnings(as.numeric(r_wb))) %>%
  filter(!is.na(r_wb))


# ---------------------------------------------------------------------------
# OECD: long-term interest rates, annual, percent per annum
oecd = read.csv(paste0("https://sdmx.oecd.org/public/rest/data/",
                       "OECD.SDD.STES,DSD_STES@DF_FINMARK,4.0/",
                       ".A.IRLT.PA......?format=csvfile")) %>%
  transmute(iso = REF_AREA, year = as.integer(TIME_PERIOD),
            r_oecd = as.numeric(OBS_VALUE)) %>% filter(!is.na(r_oecd))


# ---------------------------------------------------------------------------
# IMF: hand-downloaded extracts, found by pattern (their names carry a stamp)
imf_file = list.files("data", "IMF.STA_MFS_IR", full.names = TRUE)[1]
weo_file = list.files("data", "IMF.RES_WEO",    full.names = TRUE)[1]
if (is.na(imf_file) || is.na(weo_file))
  stop("IMF extracts missing from data/. See the README.")

# Treasury bill yields, the last fallback. Joined on name, so a few will miss.
imf_short = read.csv(imf_file, check.names = FALSE, colClasses = "character") %>%
  select(name = COUNTRY, matches("^(19|20)[0-9]{2}$")) %>%
  pivot_longer(-name, names_to = "year", values_to = "r_imf") %>%
  mutate(year = as.integer(year), r_imf = suppressWarnings(as.numeric(r_imf))) %>%
  filter(!is.na(r_imf))

# External debt is only published for the six regional aggregates
weo_ind = c(debt    = "External debt, US dollar",
            exports = "Exports of goods and services, US dollar",
            imports = "Imports of goods and services, US dollar",
            paid    = "External debt: total debt service, interest, US dollar",
            gdp     = "Gross domestic product (GDP), Current prices, US dollar")

weo_regions = read.csv(weo_file, check.names = FALSE, colClasses = "character") %>%
  filter(INDICATOR %in% weo_ind) %>%
  select(region = COUNTRY, INDICATOR, matches("^(19|20)[0-9]{2}$")) %>%
  pivot_longer(-c(region, INDICATOR), names_to = "year", values_to = "value") %>%
  mutate(year  = as.integer(year),
         value = suppressWarnings(as.numeric(value)),
         INDICATOR = names(weo_ind)[match(INDICATOR, weo_ind)]) %>%
  pivot_wider(names_from = INDICATOR, values_from = value) %>%
  transmute(region, year,
            interest = 100 * paid / debt,
            tb       = 100 * (exports - imports) / gdp,
            debt     = 100 * debt / gdp) %>%
  filter(!is.na(tb), !is.na(debt), year <= 2024) %>%
  arrange(region, year)

# One interest rate source per country, best available.
debt_interest = wb %>% inner_join(meta, by = "iso") %>%
  left_join(r_wb, by = c("iso", "year")) %>%
  left_join(oecd, by = c("iso", "year")) %>%
  left_join(imf_short, by = c("name", "year")) %>%
  mutate(r_wb = na_if(r_wb, 0),              # a zero here means missing
         tb   = 100 * tb / gdp_usd,          # as % of GDP
         debt = 100 * debt / gdp_usd) %>%
  filter(!is.na(tb), !is.na(debt)) %>%       # drop before grouping, so no group
  group_by(iso, name, region, income) %>%    # ends up empty
  arrange(year, .by_group = TRUE) %>%
  mutate(interest = if (any(!is.na(r_oecd))) r_oecd
                    else if (any(!is.na(r_wb))) r_wb
                    else r_imf,
         # most recent unbroken run of years
         block = cumsum(year != lag(year, default = first(year) - 1) + 1)) %>%
  filter(block == max(block)) %>% ungroup() %>%
  select(iso, name, region, income, year, tb, debt, interest)

write.csv(debt_interest, "data/debt_panel.csv", row.names = FALSE)
write.csv(weo_regions,   "data/weo_regions.csv",   row.names = FALSE)


# ---------------------------------------------------------------------------
# Maddison Project (the long GDP series)
mpd_raw = read.csv(paste0("https://ourworldindata.org/grapher/",
                          "gdp-per-capita-maddison-project-database.csv"))
names(mpd_raw)[4] = "gdp_pc"    # OWID calls it "GDP per capita"

maddison = mpd_raw %>%
  filter(Code %in% c("BRA", "USA"), !is.na(gdp_pc)) %>%
  select(iso = Code, year = Year, gdp_pc) %>%
  arrange(iso, year)

write.csv(maddison, "data/maddison_gdp_pc.csv", row.names = FALSE)


