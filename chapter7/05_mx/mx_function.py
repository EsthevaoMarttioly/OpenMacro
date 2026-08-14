"""
Camada 05_mx — mx_function: resolve o MX model e devolve os momentos de 2ª ordem.

Tradução de mx_function.m (Uribe). Dado (PHIM, PHIX, PSSI, RHO_TOT, STD_TOT),
resolve o sistema log-linear (gx_hx) e computa as matrizes de variância-covariância
condicionais a choques de tot (mom). Reaproveita o solver e as primitives já
construídos e validados.
"""

from pathlib import Path
import sys
import numpy as np

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG / "01_primitives"))
sys.path.insert(0, str(PKG / "04_dsge"))
from mom import mom                        # noqa: E402
from gx_hx import gx_hx, residual          # noqa: E402
from mx_model import num_eval, IDX_STATE, IDX_CONTROL   # noqa: E402


def solve(PHIM, PHIX, PSSI, RHO_TOT, STD_TOT):
    fx, fxp, fy, fyp, eta, _ = num_eval(PHIM, PHIX, PSSI, RHO_TOT, STD_TOT)
    gx, hx, flag = gx_hx(fy, fx, fyp, fxp)
    return gx, hx, eta, flag


def moments(PHIM, PHIX, PSSI, RHO_TOT, STD_TOT):
    """Variâncias condicionais a choques de tot: sigy0 (controles), sigx0 (estados)."""
    gx, hx, eta, flag = solve(PHIM, PHIX, PSSI, RHO_TOT, STD_TOT)
    varshock = eta @ eta.T
    sigy0, sigx0 = mom(gx, hx, varshock)
    return sigy0, sigx0, flag


if __name__ == "__main__":
    # calibração mediana
    gx, hx, eta, flag = solve(1.82, 1.56, 0.18, 0.5, 0.1)
    fx, fxp, fy, fyp, _, _ = num_eval(1.82, 1.56, 0.18, 0.5, 0.1)
    R = residual(gx, hx, fy, fx, fyp, fxp)
    print("exitflag (1=Blanchard-Kahn ok):", flag)
    print("nº estados:", hx.shape[0], "| resíduo max|f-cond| =", np.max(np.abs(R)))
    print("autovalores |hx| (todos <1):", np.round(np.sort(np.abs(np.linalg.eigvals(hx))), 4))

    sigy0, sigx0, _ = moments(1.82, 1.56, 0.18, 0.5, 0.1)
    # fatias de variância (variância / variância) — aqui só variâncias condicionais
    for nm in ["output_constant_prices", "c_constant_prices",
               "ivv_constant_prices", "tby_constant_prices"]:
        i = IDX_CONTROL[nm]
        print(f"  var[{nm:24s}] cond. a tot = {sigy0[i, i]:.6e}")
    print(f"  var[tot(state)] = {sigx0[IDX_STATE['tot'], IDX_STATE['tot']]:.6e}")
