# Open Macroeconomics: Mini-tasks and Exam

Authors: **Celso Nozema, Esthevão Marttioly, Mario Filho, and Victor Lucas**

Program: MSc Economics - FGV EESP

Professor: Carlos Eduardo

This repository contains the code and results for the **Open Macro** final course from FGV EESP.

All code is written with reproducibility and defensive programming in mind.

## Project Structure

```
├── code/
│   ├── main.R                           # Mini-tasks script
├── data/
│   ├── IMF_WEO_dataset.csv              # IMF-WEO Data
│   └── data.csv                         # Not yet
├── output/
├── project.RProj                        # R Project for downloading data
├── requirements.txt                     # pip install -r requirements.txt
└── README.md
```

## Data

Data was obtained in:

* [IMF: WEO Dataset](https://data.imf.org/en/datasets/IMF.RES:WEO).
** Change the name to "IMF_WEO_dataset.csv".


## Computational Environment

The analysis was conducted using R version 4.6.0 (2026-08-01) on a Windows 11 system.

## Running the project

To reproduce the analysis:

* Open the project's folder as a project.
* Open the file: code/main.py.
* Run "renv::restore()" and answer 'Y' (just the first time).
* Run the script: code/main.py.
