"""
Camada dsge — Figura 7.3: IRF da balança comercial a um choque de tot no SOE-RBC.

Tradução de edeir_ir.m (Uribe). Resolve o modelo para rho = 0, 0.25, 0.5 e mostra
a resposta da balança comercial (tb) a um choque de +1 no tot. Ilustra o efeito
ORS reforçado pelo investimento: quanto mais persistente o tot, MENOR o impacto
na tb.
"""

from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG / "01_primitives"))
from ir import ir                              # noqa: E402
from edeir_model import num_eval               # noqa: E402
from gx_hx import gx_hx                         # noqa: E402

NA = 3        # posição de a(=tot) no vetor de estados [d, r, k, a]
NTB = 6       # posição de tb no vetor de controles (0-indexado)
T = 11


def irf_tb(rho):
    fx, fxp, fy, fyp, eta, _ = num_eval(rho)
    gx, hx, flag = gx_hx(fy, fx, fyp, fxp)
    gx = np.vstack([gx, np.zeros((1, gx.shape[1]))])   # linha extra p/ tot
    gx[-1, NA] = 1.0                                    # tot = estado a
    x0 = np.zeros(hx.shape[0]); x0[NA] = 1.0            # choque de +1 no tot
    _, IRy, _ = ir(gx, hx, x0, T)                       # gx mapeia estados -> controles
    return IRy[:, NTB]


if __name__ == "__main__":
    t = np.arange(T)
    styles = {0.0: "-", 0.25: "--", 0.5: "-x"}
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for rho in (0.0, 0.25, 0.5):
        y = irf_tb(rho)
        ax.plot(t, y, styles[rho], lw=2, label=fr"$\rho={rho}$")
        print(f"rho={rho}: tb impacto={y[0]:.4f}, ano1={y[1]:.4f}")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("periods after the shock")
    ax.set_ylabel("% dev. from GDP trend")
    ax.set_title("Fig. 7.3 — IRF da balança comercial a um choque de tot (SOE-RBC)")
    ax.legend()
    fig.tight_layout()
    out = PKG / "04_dsge" / "results" / "fig_7_3_edeir_ir.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("->", out.relative_to(PKG))
