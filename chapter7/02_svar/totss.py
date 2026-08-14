"""
Camada 2 — totss: estimação do SVAR em PAINEL (pooled).

Tradução de totss.m (Uribe). Empilha os dados de todos os países num único OLS
(a defasagem é montada país-a-país com lagg ANTES de empilhar, para não misturar
o último ano de um país com o primeiro do seguinte). Usado como robustez na
Tabela 7.2 (linha "Panel Estimation").

Vetor de 5 variáveis, ordem [tot, tb, gdp, c, i] (mesma do tot_cbcs sizevar=5).

Roda para cada método de detrend:
  quadratic    -> replicação  -> 02_svar/results/
  hp, hamilton -> extensão     -> extension/svar/
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG / "01_primitives"))
from lagg import lagg                       # noqa: E402

from tot_cbcs import fit_svar, DATA, _results_dir   # noqa: E402

COLS = ["tot", "tb", "gdp", "c", "i"]


def run_panel(method, df=None):
    """Estimação em painel para um método. Devolve dict com hx, PI, R2, ssigx, v_share."""
    if df is None:
        df = pd.read_csv(DATA)
    sub = df[df.detrending_method == method]
    nv = len(COLS)

    Yb, wb = [], []
    ncou = 0
    for c in dict.fromkeys(sub.country):
        d = sub[sub.country == c].sort_values("year")[COLS].dropna()
        if len(d) < nv + 3:
            continue
        dl = lagg(d.values, 1)              # (T-1) x 2nv, por país
        Yb.append(dl[:, :nv])
        wb.append(dl[:, nv:])
        ncou += 1

    Y = np.vstack(Yb)                        # empilha contemporâneos
    w = np.vstack(wb)                        # empilha defasados
    res = fit_svar(Y, w)
    res["ncou"] = ncou
    res["nobs"] = Y.shape[0]
    res["cols"] = np.array(COLS, dtype=object)
    return res


if __name__ == "__main__":
    df = pd.read_csv(DATA)
    for method in ["quadratic", "hp", "hamilton"]:
        outdir = _results_dir(method)
        outdir.mkdir(parents=True, exist_ok=True)
        res = run_panel(method, df)
        fp = outdir / f"totss_{method}.npz"
        np.savez(fp, **res)
        vs = np.round(res["v_share"], 1)
        print(f"{method:10s}: {res['ncou']} países, {res['nobs']} obs | "
              f"v_share[tot,tb,gdp,c,i]={vs} -> {fp.relative_to(PKG)}")
