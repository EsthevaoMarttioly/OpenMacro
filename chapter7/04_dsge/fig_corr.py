"""
Camada dsge — Figura 7.4: correlação condicional entre tb e tot vs persistência.

Tradução de tot_corr_tb_tot_rho.m (Uribe). Compara:
  - linha: corr(tot, tb) condicional a choques de tot IMPLICADA PELO SOE-RBC,
    como função de rho (mostra o efeito ORS: relação negativa).
  - círculos: a mesma correlação estimada no SVAR de 2 variáveis, um por país
    (rho = hx[0,0]; corr = ssigx[0,1]/sqrt(ssigx[0,0]*ssigx[1,1])).

A conclusão do livro: o modelo prevê relação negativa, mas os dados não mostram
tal padrão -> os dados não apoiam o efeito ORS.

Os círculos dependem do método de detrend:
  quadratic -> 04_dsge/results/ ; hp, hamilton -> extension/figures/
(a linha do modelo é a mesma, independe do detrend).
"""

from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG / "01_primitives"))
from mom import mom                            # noqa: E402
from edeir_model import num_eval               # noqa: E402
from gx_hx import gx_hx                         # noqa: E402

NA = 3        # a(=tot) nos estados
NTB = 6       # tb nos controles


def model_curve(rhos):
    corr = []
    for rho in rhos:
        fx, fxp, fy, fyp, eta, _ = num_eval(rho)
        gx, hx, flag = gx_hx(fy, fx, fyp, fxp)
        if flag != 1:
            corr.append(np.nan); continue
        gx = np.vstack([gx, np.zeros((1, gx.shape[1]))])
        gx[-1, NA] = 1.0                       # tot como último controle
        nac = gx.shape[0] - 1
        varshock = eta @ eta.T
        sigy, _ = mom(gx, hx, varshock)
        cc = sigy[NTB, nac] / np.sqrt(sigy[NTB, NTB] * sigy[nac, nac])
        corr.append(cc)
    return np.array(corr)


def data_points(method):
    d = PKG / ("02_svar/results" if method == "quadratic" else "extension/svar")
    z = np.load(d / f"tot_cbcs_{method}_2.npz", allow_pickle=True)
    hx, ss = z["hx"], z["ssigx"]
    rho = hx[:, 0, 0]
    corr = ss[:, 0, 1] / np.sqrt(ss[:, 0, 0] * ss[:, 1, 1])
    return rho, corr


if __name__ == "__main__":
    rhos = np.arange(-0.5, 0.995, 0.01)
    line = model_curve(rhos)

    for method in ["quadratic", "hp", "hamilton"]:
        rho_d, corr_d = data_points(method)
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(rho_d, corr_d, "o", mfc="none", label="SVAR")
        ax.plot(rhos, line, "-", lw=2, color="tab:red", label="SOE-RBC model")
        ax.set_xlabel(r"$\rho$")
        ax.set_ylabel(r"corr$(tot, tb)$")
        ax.set_title(f"Fig. 7.4 — corr(tot,tb) condicional vs persistência — {method}")
        ax.legend(loc="lower left")
        fig.tight_layout()
        outdir = PKG / ("04_dsge/results" if method == "quadratic"
                        else "extension/figures")
        outdir.mkdir(parents=True, exist_ok=True)
        fp = outdir / f"fig_7_4_corr_{method}.png"
        fig.savefig(fp, dpi=130, bbox_inches="tight")
        plt.close(fig)
        print(f"{method:10s}: modelo corr@rho=0.5={line[np.argmin(abs(rhos-0.5))]:.3f} "
              f"-> {fp.relative_to(PKG)}")
