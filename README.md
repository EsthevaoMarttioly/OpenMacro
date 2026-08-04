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
├── output/              # Tables and figures
├── rstudio_proj.Rproj   # Open this to RStudio
├── matlab_proj.prj      # Open this to MATLAB
└── README.md
```

## Data

`data.R` is the only file that touches the network.
It pulls everything through the official R clients: `WDI` for the World Bank,
`sidrar` for IBGE/SIDRA and `rbcb` for the Banco Central's SGS.
And it writes tidy csv files into `data/`.
`main.R` reads those and never downloads anything.

| Series | Coverage | Source |
|---|---|---|
| Brazil GDP, consumption, population, tb/y | 1960+ | World Bank WDI (republished IBGE), constant LCU |
| Brazil quarterly GDP, NSA | 1996+ | [SIDRA/IBGE](https://apisidra.ibge.gov.br) table 1620, variable 583, chained volume index (1995 = 100) |
| Brazil quarterly GDP, seasonally adjusted | 1996+ | [BCB SGS](https://api.bcb.gov.br) series 22109 |
| Brazil and US real GDP per capita | 1800+ | [Maddison Project 2023](https://www.rug.nl/ggdc/historicaldevelopment/maddison/releases/maddison-project-database-2023), via Our World in Data |
| Exchange rate and openness, ~200 countries | 1960+ | World Bank WDI |

Two choices worth recording. Brazil's annual series come through WDI rather
than SIDRA because SIDRA's API only serves the quarterly system, which starts
in 1996: the pre-1996 IBGE annual data is only available through the World
Bank. And the trade balance is `tby` from the same source rather than BCB's SGS 22707,
which is the BPM6 vintage and only starts in 1995.

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
