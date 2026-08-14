"""
Camada 2 — tot_cbcs: estimação do SVAR país-a-país + decomposição de variância.

Tradução de tot_cbcs.m (Uribe). Para cada país estima, por OLS equação-a-equação,
o sistema SVAR e monta as matrizes hx (dinâmica) e PI (impacto dos choques),
além de R2, variâncias condicionais a choques de tot (ssigx) e a participação da
variância explicada por tot (v_share).

sizevar = 2  -> SVAR de 2 variáveis (tot, tb)          [seção 7.1-7.4]
sizevar = 5  -> SVAR de 5 variáveis (tot, tb, gdp, c, i) [seção "How Important..."]

Roda para cada método de detrend do dataset (quadratic, hp, hamilton).
  quadratic -> replicação do livro  -> resultados em 02_svar/results/
  hp, hamilton -> extensão           -> resultados em extension/svar/

Sistema estimado (ver slides): com A0 x_t = A x_{t-1} + PI1 eps,
  hx = A0^{-1} A ,  PI = A0^{-1} PI1 ,  PI1 = chol(cov(u), 'lower').
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG / "01_primitives"))
from lagg import lagg                                  # noqa: E402
from mom import mom                                    # noqa: E402
from variance_decomposition import variance_decomposition  # noqa: E402

DATA = PKG / "00_data" / "replication_detrended.csv"
COLS = {2: ["tot", "tb"], 5: ["tot", "tb", "gdp", "c", "i"]}


def fit_svar(Y, w):
    """Núcleo da estimação SVAR a partir dos blocos contemporâneo (Y) e defasado (w).

    Y, w : (N x nv). Usado tanto país-a-país (tot_cbcs) quanto em painel (totss),
    bastando empilhar os blocos de vários países antes de chamar.
    Devolve dict com hx, PI, R2, sigx, ssigx, v_share.
    """
    Y = np.asarray(Y, float)
    w = np.asarray(w, float)
    nv = Y.shape[1]
    X = np.hstack([w, np.ones((w.shape[0], 1))])   # [defasados, constante]

    u = np.zeros((X.shape[0], nv))
    A = np.zeros((nv, nv))
    A0 = np.eye(nv)

    # --- equação 1: tot_t sobre [tot_{t-1}, const] ---
    y = Y[:, 0]
    w1 = np.column_stack([X[:, 0], X[:, -1]])
    aa, *_ = np.linalg.lstsq(w1, y, rcond=None)
    A[0, 0] = aa[0]
    u[:, 0] = y - w1 @ aa

    # --- equações 2..nv: cada var sobre [tot_t, defasados, const] ---
    if nv > 1:
        Xy = np.column_stack([y, X])            # [tot_t, defasados(nv), const]
        aa, *_ = np.linalg.lstsq(Xy, Y[:, 1:], rcond=None)
        A[1:, :] = aa[1:nv + 1, :].T            # coefs dos defasados -> A
        A0[1:, 0] = -aa[0, :].T                 # coef contemporâneo de tot -> A0
        u[:, 1:] = Y[:, 1:] - Xy @ aa

    R2 = 1.0 - u.var(axis=0, ddof=1) / Y.var(axis=0, ddof=1)
    SIGMA = np.cov(u, rowvar=False)             # normaliza por N-1 (como MATLAB)
    PI1 = np.linalg.cholesky(np.atleast_2d(SIGMA))   # lower
    PI = np.linalg.solve(A0, PI1)
    hx = np.linalg.solve(A0, A)

    _, sigx = mom(np.eye(nv), hx, PI @ PI.T)         # variância total
    PIz = PI.copy(); PIz[:, 1:] = 0.0                # só choque de tot
    _, ssigx = mom(np.eye(nv), hx, PIz @ PIz.T)      # variância condicional
    _, vdx, _, _ = variance_decomposition(np.eye(nv), hx, PI)
    v_share = vdx[0, :] * 100.0                       # % explicado por tot

    return dict(hx=hx, PI=PI, R2=R2, sigx=sigx, ssigx=ssigx, v_share=v_share)


def estimate_country(D):
    """Estima o SVAR de um país. D: (T x nv), colunas na ordem [tot, tb, ...]."""
    D = np.asarray(D, float)
    nv = D.shape[1]
    dl = lagg(D, 1)                 # (T-1) x 2nv: [contemporâneo | defasado]
    return fit_svar(dl[:, :nv], dl[:, nv:])


def run(method, sizevar, df=None):
    """Estima todos os países para um (método, sizevar). Devolve arrays empilhados."""
    if df is None:
        df = pd.read_csv(DATA)
    cols = COLS[sizevar]
    sub = df[df.detrending_method == method]
    countries = list(dict.fromkeys(sub.country))       # ordem preservada
    nv = len(cols)
    res = {k: [] for k in ["hx", "PI", "R2", "sigx", "ssigx", "v_share"]}
    used = []
    for c in countries:
        d = sub[sub.country == c].sort_values("year")[cols]
        d = d.dropna()                                  # hamilton perde anos iniciais
        if len(d) < nv + 3:                             # amostra mínima p/ estimar
            continue
        r = estimate_country(d.values)
        for k in res:
            res[k].append(r[k])
        used.append(c)
    out = {k: np.array(v) for k, v in res.items()}      # ex.: hx -> (ncou, nv, nv)
    out["countries"] = np.array(used, dtype=object)
    out["cols"] = np.array(cols, dtype=object)
    return out


def _results_dir(method):
    if method == "quadratic":
        return PKG / "02_svar" / "results"
    return PKG / "extension" / "svar"


if __name__ == "__main__":
    df = pd.read_csv(DATA)
    for method in ["quadratic", "hp", "hamilton"]:
        outdir = _results_dir(method)
        outdir.mkdir(parents=True, exist_ok=True)
        for sizevar in (2, 5):
            out = run(method, sizevar, df)
            fp = outdir / f"tot_cbcs_{method}_{sizevar}.npz"
            np.savez(fp, **out)
            print(f"{method:10s} sizevar={sizevar}: {len(out['countries'])} países "
                  f"-> {fp.relative_to(PKG)}")
