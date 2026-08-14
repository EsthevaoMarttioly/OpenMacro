"""
Camada 05_mx — Modelo simbólico MX (dois setores) e derivadas.

Tradução de mx_model.m (Uribe) com sympy. 34 equações E_t f(yp,y,xp,x)=0.
Armington Cobb-Douglas (ces=0, MU=1). Aplica a regra log (d/d log v = v·d/dv) e
avalia no estado estacionário.

Estados   x = [d, km, kx, tfp_m, tfp_x, tot]           (todos em log menos d)
Controles y = [c, hm, h_x, im, ix, am, ax, px, pm, ym, yx, wm, wx, ux, um,
               m, x, output, la, kmfu, kxfu, ivv,
               output_constant_prices, c_constant_prices, ivv_constant_prices,
               tby, r, tby_constant_prices]            (os 3 últimos em nível)

O estado estacionário independe de PHIM/PHIX/PSSI/RHO_TOT; logo calculamos o ss
uma vez e deixamos esses 4 parâmetros livres nas derivadas (lambdify).
"""

import numpy as np
import sympy as sp

from mx_ss import steady_state

# ---- nomes ----
STATES = ["d", "km", "kx", "tfp_m", "tfp_x", "tot"]
STATE_LOG = {"d": False, "km": True, "kx": True, "tfp_m": True, "tfp_x": True, "tot": True}
CONTROLS = ["c", "hm", "h_x", "im", "ix", "am", "ax", "px", "pm", "ym", "yx",
            "wm", "wx", "ux", "um", "m", "x", "output", "la", "kmfu", "kxfu",
            "ivv", "output_constant_prices", "c_constant_prices",
            "ivv_constant_prices", "tby", "r", "tby_constant_prices"]
CONTROL_LOG = {n: True for n in CONTROLS}
for n in ["tby", "r", "tby_constant_prices"]:
    CONTROL_LOG[n] = False
NAMES = STATES + CONTROLS

cur = {n: sp.Symbol(n) for n in NAMES}
nxt = {n: sp.Symbol(n + "p") for n in NAMES}

# parâmetros livres nas derivadas
PHIM, PHIX, PSSI, RHO_TOT = sp.symbols("PHIM PHIX PSSI RHO_TOT")
# parâmetros fixos (símbolos; substituídos numericamente depois)
SIGG, OMEGAM, OMEGAX, ALFAM, ALFAX, BETTA, DELTA, CHI, MU, RSTAR, DBAR, \
    TOT, TFP_X, TFP_M, OUTPUT, PX, PM = sp.symbols(
        "SIGG OMEGAM OMEGAX ALFAM ALFAX BETTA DELTA CHI MU RSTAR DBAR "
        "TOT TFP_X TFP_M OUTPUT PX PM")

# atalhos de variáveis
V = lambda n: cur[n]
Vp = lambda n: nxt[n]
(c, hm, h_x, im, ix, am, ax, px, pm, ym, yx, wm, wx, ux, um, m, x, output,
 la, kmfu, kxfu, ivv, ocp, ccp, ivcp, tby, r, tbycp) = [V(n) for n in CONTROLS]
d, km, kx, tfp_m, tfp_x, tot = [V(n) for n in STATES]
(cp, hmp, h_xp, imp, ixp, amp, axp, pxp, pmp, ymp, yxp, wmp, wxp, uxp, ump,
 mp, xp, outputp, lap, kmfup, kxfup, ivvp, ocpp, ccpp, ivcpp, tbyp, rp,
 tbycpp) = [Vp(n) for n in CONTROLS]
dp, kmp, kxp, tfp_mp, tfp_xp, totp = [Vp(n) for n in STATES]

# ---- utilidade e agregador ----
U = ((c - hm**OMEGAM/OMEGAM - h_x**OMEGAX/OMEGAX)**(1 - SIGG) - 1) / (1 - SIGG)
A = am**CHI * ax**(1 - CHI)                        # Cobb-Douglas (ces=0)
Uc, Uhm, Uhx = U.diff(c), U.diff(hm), U.diff(h_x)

# ---- 34 equações ----
f = sp.Matrix([
    -la + Uc,                                                          # e1
    wm + Uhm / Uc,                                                     # e2
    wx + Uhx / Uc,                                                     # e3
    -la + BETTA*(1 + r)*lap,                                           # e4
    -la*(1 + PHIM*(kmp - km)) + BETTA*lap*(ump + 1 - DELTA + PHIM*(kmfup - kmp)),  # e5
    -la*(1 + PHIX*(kxp - kx)) + BETTA*lap*(uxp + 1 - DELTA + PHIX*(kxfup - kxp)),  # e6
    -kmp + (1 - DELTA)*km + im,                                        # e7
    -kxp + (1 - DELTA)*kx + ix,                                        # e8
    -pm + A.diff(am),                                                  # e9
    -px + A.diff(ax),                                                  # e10
    -tfp_m*km**ALFAM*hm**(1 - ALFAM) + ym,                             # e11
    -tfp_x*kx**ALFAX*h_x**(1 - ALFAX) + yx,                            # e12
    -um + pm*tfp_m*ALFAM*km**(ALFAM - 1)*hm**(1 - ALFAM),              # e13
    -ux + px*tfp_x*ALFAX*kx**(ALFAX - 1)*h_x**(1 - ALFAX),             # e14
    -wm + pm*tfp_m*(1 - ALFAM)*km**ALFAM*hm**(-ALFAM),                 # e15
    -wx + px*tfp_x*(1 - ALFAX)*kx**ALFAX*h_x**(-ALFAX),                # e16
    -A + c + im + ix + PHIX/2*(kxp - kx)**2 + PHIM/2*(kmp - km)**2,    # e17
    pm*(-am + ym) + m,                                                 # e18
    px*(-ax + yx) - x,                                                 # e19
    -d - m + dp/(1 + r) + x,                                           # e20
    -r + RSTAR + PSSI*(sp.exp(dp - DBAR) - 1),                         # e21
    -tot + px/pm,                                                      # e22
    -sp.log(totp/TOT) + RHO_TOT*sp.log(tot/TOT),                       # e23a
    -sp.log(tfp_xp/TFP_X) + RHO_TOT*sp.log(tfp_x/TFP_X),               # e23b
    -sp.log(tfp_mp/TFP_M) + RHO_TOT*sp.log(tfp_m/TFP_M),               # e23c
    -kxfu + kxp,                                                       # e24a
    -kmfu + kmp,                                                       # e24b
    -output + px*yx + pm*ym,                                           # e25a
    -ivv + im + ix,                                                    # e25b
    -tby + (x - m)/OUTPUT,                                             # e26
    -ocp + PX*yx + PM*ym,                                              # e27a
    -ccp + c/output*ocp,                                               # e27b
    -ivcp + ivv/output*ocp,                                            # e27c
    -tbycp + (x - m)/output*ocp/OUTPUT,                                # e27d
])


def _deriv(names, logflags, isp):
    cols = []
    for n in names:
        v = nxt[n] if isp else cur[n]
        col = f.diff(v)
        if logflags[n]:
            col = col * v
        cols.append(col)
    return sp.Matrix.hstack(*cols)


_FX = _deriv(STATES, STATE_LOG, False)
_FXP = _deriv(STATES, STATE_LOG, True)
_FY = _deriv(CONTROLS, CONTROL_LOG, False)
_FYP = _deriv(CONTROLS, CONTROL_LOG, True)


def _ss_subs():
    ss, consts = steady_state()          # ss independe de PHIM/PHIX/PSSI/RHO_TOT/STD_TOT
    subs = {}
    for n in NAMES:
        subs[cur[n]] = ss[n]
        subs[nxt[n]] = ss[n]
    for pk in ["SIGG", "OMEGAM", "OMEGAX", "ALFAM", "ALFAX", "BETTA", "DELTA",
               "CHI", "MU", "RSTAR", "DBAR", "TOT", "TFP_X", "TFP_M", "OUTPUT",
               "PX", "PM"]:
        subs[sp.Symbol(pk)] = consts[pk]
    return subs


_SUB = _ss_subs()
_args = (PHIM, PHIX, PSSI, RHO_TOT)
_fx_l = sp.lambdify(_args, _FX.subs(_SUB), "numpy")
_fxp_l = sp.lambdify(_args, _FXP.subs(_SUB), "numpy")
_fy_l = sp.lambdify(_args, _FY.subs(_SUB), "numpy")
_fyp_l = sp.lambdify(_args, _FYP.subs(_SUB), "numpy")
_f_l = sp.lambdify(_args, f.subs(_SUB), "numpy")

IDX_STATE = {n: i for i, n in enumerate(STATES)}
IDX_CONTROL = {n: i for i, n in enumerate(CONTROLS)}


def num_eval(PHIM_v, PHIX_v, PSSI_v, RHO_TOT_v, STD_TOT_v):
    a = (PHIM_v, PHIX_v, PSSI_v, RHO_TOT_v)
    fx = np.array(_fx_l(*a), float)
    fxp = np.array(_fxp_l(*a), float)
    fy = np.array(_fy_l(*a), float)
    fyp = np.array(_fyp_l(*a), float)
    eta = np.zeros((len(STATES), 1)); eta[IDX_STATE["tot"], 0] = STD_TOT_v
    fss = np.array(_f_l(*a), float).ravel()
    return fx, fxp, fy, fyp, eta, fss


if __name__ == "__main__":
    fx, fxp, fy, fyp, eta, fss = num_eval(1.82, 1.56, 0.18, 0.5, 0.1)
    print("dims: fx", fx.shape, "fy", fy.shape, "fxp", fxp.shape, "fyp", fyp.shape)
    print("nº equações:", f.shape[0], "| estados:", len(STATES), "| controles:", len(CONTROLS))
    print("resíduo estado estacionário max|f(ss)| =", np.max(np.abs(fss)))
