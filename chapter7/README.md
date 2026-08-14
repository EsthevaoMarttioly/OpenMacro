# Terms of Trade — replicação em Python

Replicação da parte empírica do Capítulo 7 ("Importable Goods, Exportable Goods,
and the Terms of Trade") de Uribe & Schmitt-Grohé, *Open Economy Macroeconomics*,
Princeton University Press, 2017.

## Estrutura

| Pasta          | Papel                              | Espelha (MATLAB)                                   |
|----------------|------------------------------------|---------------------------------------------------|
| `data/`        | dados base + conversão do `.mat`   | `usg_tot_data.mat`                                |
| `data/processed/` | `.mat`/`.npz` intermediários     | `tot_cbcs2.mat`, `tot_cbcs5.mat`                  |
| `primitives/`  | funções utilitárias (lego)         | `lagg.m`, `ir.m`, `mom.m`, `variance_decomposition.m` |
| `svar/`        | motor econométrico                 | `tot_cbcs.m`, `totss.m`                           |
| `dsge/`        | modelo SOE-RBC + solver 1ª ordem   | `edeir_*.m`, `gx_hx`                              |
| `figures/`     | scripts que geram as figuras       | `tot_ir_cbc.m`, `ors_test.m`, `edeir_ir.m`, ...  |
| `tables/`      | scripts que geram as tabelas       | `table_ar1_cbc.m`, `table_v_share_cbc.m`         |
| `outputs/`     | resultados gerados (.png/.csv)     | (saída de tela/console no MATLAB)                |

## Ordem de execução

1. `data/` — converter `usg_tot_data.mat` para formato Python
2. `svar/tot_cbcs.py` — rodar com sizevar=2 e sizevar=5 (gera os intermediários)
3. `tables/` e `figures/` — consomem os intermediários

`primitives/` é chamado automaticamente pelas camadas acima.

## Mapa figuras/tabelas

- Tabela 7.1 → `tables/table_7_1_ar1.py`
- Tabela 7.2 → `tables/table_7_2_vshare.py`
- Figura 7.1 / 7.5 → `figures/fig_7_1_5_ir.py`
- Figura 7.2 → `figures/fig_7_2_ors.py`
- Figura 7.3 → `figures/fig_7_3_edeir_ir.py`
- Figura 7.4 → `figures/fig_7_4_corr.py`
