"""
Camada 05_mx — Figura de IRFs do MX model.

Tradução de mx_ir.m + plot_mx_ir.m (Uribe). Para cada país resolve o modelo com a
calibração de phi_pssi_cbc.mat, calcula a resposta a um choque de +10% no tot, e
plota a MEDIANA entre países num painel de 18 variáveis.
"""

from pathlib import Path
import sys
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG / "01_primitives"))
from ir import ir                                  # noqa: E402
from mx_function import solve                       # noqa: E402
from mx_model import IDX_STATE, IDX_CONTROL, STATES, CONTROLS   # noqa: E402
from mx_ss import steady_state                      # noqa: E402

MX = Path(__file__).resolve().parents[2] / "usg_mx"
TIR = 11

# (título, tipo, nome)  tipo: 'c'=controle, 's'=estado, 'd'=dívida/OUTPUT
PANELS = [
    ("Terms of Trade", "s", "tot"), ("Imports", "c", "m"), ("Exports", "c", "x"),
    ("Trade Balance", "c", "tby_constant_prices"),
    ("Output in Import Sector", "c", "ym"), ("Output in Export Sector", "c", "yx"),
    ("Investment", "c", "ivv_constant_prices"),
    ("Investment in Import Sector", "c", "im"),
    ("Investment in Export Sector", "c", "ix"),
    ("GDP", "c", "output_constant_prices"),
    ("Employment in Import Sector", "c", "hm"),
    ("Employment in Export Sector", "c", "h_x"),
    ("Consumption", "c", "c_constant_prices"),
    ("Absorption of Importables", "c", "am"),
    ("Absorption of Exportables", "c", "ax"),
    ("External Debt", "d", "d"),
    ("Wage in Import Sector", "c", "wm"),
    ("Wage in Export Sector", "c", "wx"),
]


def all_irfs():
    phi = sio.loadmat(MX / "phi_pssi_cbc.mat")
    PHIM = phi["PHIM_cbc"].ravel(); PHIX = phi["PHIX_cbc"].ravel()
    PSSI = phi["PSSI_cbc"].ravel(); RHO = phi["RHO_TOT_cbc"].ravel()
    STD = phi["STD_TOT_cbc"].ravel()
    ncou = len(PHIM)
    _, consts = steady_state()
    OUTPUT = consts["OUTPUT"]

    ny, nx = len(CONTROLS), len(STATES)
    IRy = np.zeros((ncou, TIR, ny))
    IRx = np.zeros((ncou, TIR, nx))
    for k in range(ncou):
        gx, hx, eta, flag = solve(PHIM[k], PHIX[k], PSSI[k], RHO[k], STD[k])
        x0 = np.zeros(nx); x0[IDX_STATE["tot"]] = 10.0     # +10% no tot
        _, iry, irx = ir(gx, hx, x0, TIR)
        IRy[k], IRx[k] = iry, irx
    return np.median(IRy, axis=0), np.median(IRx, axis=0), OUTPUT


if __name__ == "__main__":
    IRy, IRx, OUTPUT = all_irfs()
    t = np.arange(TIR)
    fig, axes = plt.subplots(6, 3, figsize=(11, 15))
    axes = axes.ravel()
    for ax, (title, typ, name) in zip(axes, PANELS):
        if typ == "s":
            y = IRx[:, IDX_STATE[name]]
        elif typ == "d":
            y = IRx[:, IDX_STATE[name]] / OUTPUT
        else:
            y = IRy[:, IDX_CONTROL[name]]
        ax.plot(t, y, lw=2.5)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(title, fontsize=10)
    fig.suptitle("MX Model — IRFs medianas a um choque de +10% no tot", y=1.002)
    fig.tight_layout()
    out = PKG / "05_mx" / "results" / "fig_mx_ir.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print("tot impacto:", round(IRx[0, IDX_STATE['tot']], 3),
          "| tb impacto:", round(IRy[0, IDX_CONTROL['tby_constant_prices']], 4))
    print("->", out.relative_to(PKG))
