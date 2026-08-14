Replication files for regressions in  ``Country Spreads and Emerging Countries: Who Drives Whom?,' by Martín Uribe and  Vivian Z. Yue, Journal of International Economics, 69, June 2006, 6-36. 

- Replication SGU chapter 6 - SLIDES.pdf: slide set that summarizes Chapter 6 () of Open Economy Macroeconomics (2017) by Schmitt-Grohé, S., & Uribe, M.

- Open_macro_2026___final_exam___cap6_replications.pdf: same information as the slide set organized as a short paper

- uribe_yue_derivation.pdf: Guide with all the derivations to the model equations included in the Dynare code

FOLDER>> Replication Uribe Yue 2006\

- statadata.xlsx: original panel data

- uribe_yue_VAR_baseline_v2.mlx: runs the baseline SVAR with the original functional form and data, replicating results

- replicate_IRF_uribe_yue_matching.m: calls the Dynare model with the necessary inputs to run the estimation

- DSGE_uribe_yue_IRFmatch.MOD: baseline estimated model converted to Dynare´s .mod file

FOLDER>> SGU_cap6_extension_Brazil\

- Brazil_extended_data.xlsx: compatible tiem series data extended from 1996 to 2026 for Brazil, following analogous specification as the original dataset from the authors

- uribe_yue_VAR_Brazil.mlx: runs the single-country SVAR with functional form analogous to the original SVAR

- replicate_IRF_uribe_yue_matching_Brazil.m: calls the Dynare model with the necessary inputs to run the estimation

- DSGE_uribe_yue_IRFmatch_Brazil.MOD:  single-country estimated model converted to Dynare´s .mod file (identical to DSGE_uribe_yue_IRFmatch.MOD)


