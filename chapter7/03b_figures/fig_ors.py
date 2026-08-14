"""
Camada 3b — Figura 7.2: teste do efeito ORS.

Tradução de ors_test.m (Uribe). Para cada país plota a persistência do tot
(rho = hx[0,0]) contra a resposta de IMPACTO da balança comercial a um choque de
+10% no tot (10 * PI[1,0] / PI[0,0]). Se o efeito ORS estivesse nos dados,
veríamos um padrão negativo — o livro conclui que não há.

Usa o SVAR de 2 variáveis (sizevar=2).
quadratic -> 03b_figures/results/ ; hp, hamilton -> extension/figures/
"""

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PKG = Path(__file__).resolve().parents[1]


def _npz(method):
    d = PKG / ("02_svar/results" if method == "quadratic" else "extension/svar")
    return np.load(d / f"tot_cbcs_{method}_2.npz", allow_pickle=True)


def ors_points(method):
    z = _npz(method)
    rho = z["hx"][:, 0, 0]
    impact = 10.0 * z["PI"][:, 1, 0] / z["PI"][:, 0, 0]
    return rho, impact, z["countries"]


def plot(method, outdir):
    rho, impact, _ = ors_points(method)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(rho, impact, "o", mfc="none")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel(r"$\rho$")
    ax.set_ylabel("% dev. from GDP trend")
    ax.set_title(f"Teste ORS: impacto de tot na tb vs persistência — {method}")
    # linha de regressão (para quantificar a relação — útil na extensão)
    b1, b0 = np.polyfit(rho, impact, 1)
    xs = np.linspace(rho.min(), rho.max(), 50)
    ax.plot(xs, b0 + b1 * xs, "--", color="tab:red", lw=1.5,
            label=f"inclinação = {b1:.2f}")
    ax.legend()
    fig.tight_layout()
    fp = outdir / f"fig_7_2_ors_{method}.png"
    fig.savefig(fp, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return fp, b1


if __name__ == "__main__":
    for method in ["quadratic", "hp", "hamilton"]:
        outdir = PKG / ("03b_figures/results" if method == "quadratic"
                        else "extension/figures")
        outdir.mkdir(parents=True, exist_ok=True)
        fp, slope = plot(method, outdir)
        print(f"{method:10s}: inclinação da regressão ORS = {slope:+.3f} "
              f"-> {fp.relative_to(PKG)}")
