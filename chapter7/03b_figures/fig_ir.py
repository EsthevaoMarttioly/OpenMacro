"""
Camada 3b — Figuras 7.1 e 7.5: respostas a impulso a um choque de +10% no tot.

Tradução de tot_ir_cbc.m (Uribe). Para cada país constrói a IRF a partir de
(hx, PI) e do impulso x0 = PI[:,0]/PI[0,0]*10 (normalizado para tot +10% no
impacto); a figura mostra a MÉDIA entre países.

  sizevar=2 -> Figura 7.1 (tot, tb)
  sizevar=5 -> Figura 7.5 (tot, tb, output, consumo, investimento)

quadratic -> 03b_figures/results/ ; hp, hamilton -> extension/figures/
"""

from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG / "01_primitives"))
from ir import ir                                   # noqa: E402

T = 11
TITLES = {2: ["Terms of Trade", "Trade Balance"],
          5: ["Terms of Trade", "Trade Balance", "Output", "Consumption", "Investment"]}
YLAB = {0: "% dev. from trend", 1: "% dev. from GDP trend"}


def _npz(method, sizevar):
    d = PKG / ("02_svar/results" if method == "quadratic" else "extension/svar")
    return np.load(d / f"tot_cbcs_{method}_{sizevar}.npz", allow_pickle=True)


def mean_irf(method, sizevar):
    z = _npz(method, sizevar)
    hx, PI = z["hx"], z["PI"]
    nv = sizevar
    ncou = hx.shape[0]
    IRs = np.zeros((ncou, T, nv))
    for k in range(ncou):
        P, HX = PI[k], hx[k]
        x0 = P[:, 0] / P[0, 0] * 10.0        # impulso: tot +10% no impacto
        _, IRy, _ = ir(np.eye(nv), HX, x0, T)   # IRy = resposta das nv variáveis
        IRs[k] = IRy
    return IRs.mean(axis=0)                   # (T, nv)


def plot(method, sizevar, outdir):
    pIR = mean_irf(method, sizevar)
    t = np.arange(T)
    titles = TITLES[sizevar]
    if sizevar == 2:
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
        order = [0, 1]
    else:
        fig, axes = plt.subplots(3, 2, figsize=(9, 9))
        axes = axes.ravel()
        order = [0, 1, 2, 3, 4]
    for j, ax_i in enumerate(order):
        ax = axes[ax_i]
        ax.plot(t, pIR[:, j], lw=3)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(titles[j])
        ax.set_xlabel("years after the shock")
        ax.set_ylabel(YLAB.get(1 if j == 1 else 0))
    if sizevar == 5:
        axes[5].axis("off")
    fig.suptitle(f"IRF a +10% no tot — {method}", y=1.0)
    fig.tight_layout()
    fignum = "7_1" if sizevar == 2 else "7_5"
    fp = outdir / f"fig_{fignum}_ir_{method}.png"
    fig.savefig(fp, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return fp, pIR


if __name__ == "__main__":
    for method in ["quadratic", "hp", "hamilton"]:
        outdir = PKG / ("03b_figures/results" if method == "quadratic"
                        else "extension/figures")
        outdir.mkdir(parents=True, exist_ok=True)
        for sizevar in (2, 5):
            fp, pIR = plot(method, sizevar, outdir)
            print(f"{method:10s} sizevar={sizevar}: tot impacto={pIR[0,0]:.2f} "
                  f"tb impacto={pIR[0,1]:.3f} -> {fp.relative_to(PKG)}")
