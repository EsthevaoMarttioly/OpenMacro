# Open Macroeconomics: Mini-tasks and Exam

Authors: **Celso Nozema, Esthevão Marttioly, Mario Filho, and Victor Lucas**

Program: MSc Economics - FGV EESP

Professor: Carlos Eduardo

This repository contains the code and results for the **Open Macro** final course from FGV EESP.

## Project Structure

```
├── code/
│   ├── data.R           # Downloads and tidies every series
│   ├── main.R           # Q1 to Q6
│   └── results.R        # Tables and figures, called by main.R
├── dynare/              # Q7, Q8, Q9
├── data/                # Built by data.R
├── output/
│   ├── figures/         # png, from R and from Dynare
│   └── tables/          # csv and tex
├── uribe yue replication/   # Final exam, Chapter 6
├── OpenMacro_Chapter6.pdf   # Course slides
├── OpenMacro_Chapter6_Summary.pdf
├── OpenMacro_Chapter7.pdf
├── OpenMacro_Tasks.pdf      # The questions
├── rstudio_proj.Rproj   # Open this to RStudio
├── matlab_proj.prj      # Open this to MATLAB
└── README.md
```

## Final exam: Uribe and Yue (2006)

`uribe yue replication/` replicates the regressions in Uribe and Yue,
"Country Spreads and Emerging Countries: Who Drives Whom?", *Journal of
International Economics* 69, 2006, 6-36. It has its own `readme.txt`.

```
uribe yue replication/
├── Replication Uribe Yue 2006/      # Original data and specification
│   ├── statadata.xls                        # Original panel
│   ├── uribe_yue_VAR_baseline_v2.mlx        # Baseline SVAR
│   ├── replicate_IRF_uribe_yue_matching.m   # Calls Dynare for the estimation
│   └── DSGE_uribe_yue_IRFmatch.mod          # Estimated model
├── SGU_cap6_extension_Brazil/       # Same exercise, Brazil 1996 to 2026
│   ├── Brazil_extended_data.xlsx
│   ├── uribe_yue_VAR_Brazil.mlx
│   ├── replicate_IRF_uribe_yue_matching_Brazil.m
│   └── DSGE_uribe_yue_IRFmatch_Brazil.mod
├── uribe_yue_derivation.pdf         # Derivations behind the .mod equations
├── Archive/                         # Earlier model versions, v0 to v6
└── readme.txt
```

Run the `.mlx` live script first for the SVAR, then the `replicate_*.m` script,
which calls Dynare and does the IRF matching. The Brazil folder mirrors the
baseline: the `.mod` file is identical, only the data and the country differ.

## Data

`data.R` is the only file that touches the network.
It pulls everything through the official R clients: `WDI` for the World Bank,
`sidrar` for IBGE/SIDRA and `rbcb` for the Banco Central's SGS, plus the OECD
SDMX endpoint.
And it writes tidy csv files into `data/`.
`main.R` reads those and never downloads anything.

| Series | Coverage | Source |
|---|---|---|
| Brazil GDP, consumption, population, tb/y | 1960+ | World Bank WDI (republished IBGE), constant LCU |
| Brazil quarterly GDP, NSA | 1996+ | [SIDRA/IBGE](https://apisidra.ibge.gov.br) table 1620, variable 583, chained volume index (1995 = 100) |
| Brazil quarterly GDP, seasonally adjusted | 1996+ | [BCB SGS](https://api.bcb.gov.br) series 22109 |
| Brazil and US real GDP per capita | 1800+ | [Maddison Project 2023](https://www.rug.nl/ggdc/historicaldevelopment/maddison/releases/maddison-project-database-2023), via Our World in Data |
| Exchange rate and openness, ~200 countries | 1960+ | World Bank WDI |
| Trade balance, external debt stock, interest on new debt | 1970+ | World Bank WDI |
| Long-term interest rates | 1980+ | [OECD SDMX](https://sdmx.oecd.org), `DSD_STES@DF_FINMARK`, `IRLT` |
| Treasury bill yields, external debt by region | 1950+ | IMF, downloaded by hand |

Two choices worth recording. Brazil's annual series come through WDI rather
than SIDRA because SIDRA's API only serves the quarterly system, which starts
in 1996: the pre-1996 IBGE annual data is only available through the World
Bank. And the trade balance is `tby` from the same source rather than BCB's SGS 22707,
which is the BPM6 vintage and only starts in 1995.

### Data downloaded by hand

Three extracts have no usable API. Drop the csv into `data/` without renaming:
`data.R` finds them by pattern.

From [data.imf.org](https://data.imf.org):

* **WEO**, dataset `IMF.RES:WEO`. The six regional aggregates, with external
  debt, exports and imports of goods and services, external debt interest paid,
  and GDP at current prices.
* **MFS**, dataset `IMF.STA:MFS_IR`. All countries,
  "Government securities: Treasury bills yields".

From [databank.worldbank.org](https://databank.worldbank.org), International
Debt Statistics, saved as `Debt_Interest.csv`:

* `DT.DOD.DECT.CD` and `DT.INR.DPPG`, all countries, all years. `DT.INR.DPPG`
  is in IDS only, so neither the `WDI` package nor the plain API returns it.

### Interest rates

One source per country, best available, rather than splicing definitions year
by year. Long-term (10 year) rates from the OECD where they exist. For emerging
countries, mostly absent from the OECD, the average interest on new external
debt commitments from World Bank IDS. Treasury bill yields from the IMF as a
last resort. The World Bank has no external debt for many developing countries,
so the Q1 panel covers fewer of them than the country list suggests.

## Computational Environment

R 4.6.0 on Windows 11. Package versions are locked with `renv`.

## Running the project

### R, questions 1 to 6

* Open `rstudio_proj.Rproj` in RStudio. This sets the working directory to the
  repository root, which every path in the code assumes.
* Run `renv::restore()` and answer 'Y' (just the first time).
* Run `code/main.R`.

`main.R` calls `data.R` automatically if `data/` is empty.
Run `data.R` by hand to refresh the data.

### MATLAB and Dynare, questions 7 to 9

* Open `matlab_proj.prj` to open MATLAB. This opens the repository as a
  project and puts the folders on the path.
* Move to `dynare/` and run `dynare q7_capital_hump`, `dynare q8_ghh_sigma`
  and `dynare q9_rho_near_one`.

### MATLAB and Dynare, final exam

* Same project file.
* In `uribe yue replication/Replication Uribe Yue 2006/`, open
  `uribe_yue_VAR_baseline_v2.mlx` and run it, then run
  `replicate_IRF_uribe_yue_matching.m`.
* For the Brazilian extension, the same two steps in
  `SGU_cap6_extension_Brazil/` with the `_Brazil` files.

## References

Uribe, M. and Yue, V. Z. (2006). Country spreads and emerging countries: who
drives whom? *Journal of International Economics*, 69, 6-36.

Schmitt-Grohé, S. and Uribe, M. (2017). *Open Economy Macroeconomics*.
Princeton University Press. Chapter 6.
