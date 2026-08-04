#=
# ---------------------------------------------------------------------------
#  brazil_gdp_quarterly.csv     volume index, NSA                 1996+   SIDRA
#  brazil_gdp_quarterly_sa.csv  volume index, seasonally adj.     1996+   BCB
#  brazil_annual.csv            GDP, C, pop, tb/y                 1960+   WDI
#  country_panel.csv            exchange rate, openness           1960+   WDI
#  maddison_gdp_pc.csv          real GDP per capita, BRA and USA  1800+   MPD
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
           open = "NE.TRD.GNFS.ZS")  # (X + M) / Y, %

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
  select(iso, name, region, income, year, er, open) %>%
  arrange(iso, year)

write.csv(brazil, "data/brazil_annual.csv", row.names = FALSE)
write.csv(panel,  "data/country_panel.csv", row.names = FALSE)


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


