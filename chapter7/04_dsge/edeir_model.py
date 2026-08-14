"""
Camada dsge — Modelo simbólico EDEIR e derivadas (fx, fxp, fy, fyp).

Tradução de edeir_model.m (Uribe), usando sympy no lugar do Symbolic Math Toolbox.
Monta as 15 equações de equilíbrio E_t f(yp, y, xp, x) = 0, aplica a aproximação
log-linear (variáveis "em log" recebem a regra da cadeia: d/d log(v) = v * d/dv) e
avalia as derivadas no estado estacionário.

Estados   x = [d, r, k, a]                         (k, a em log)
Controles y = [c, ivv, output, h, la, kfu,         (esses 6 em log)
               tb, tby, ca, cay, tb_o_tot]         (esses 5 em nível)

Devolve as matrizes numéricas fx (15x4), fxp (15x4), fy (15x11), fyp (15x11) e
o vetor ETASHOCK (4x1). Validado contra edeir_num_eval.m (gabarito do MATLAB).
"""

import numpy as np
import sympy as sp

from edeir_ss import PARAMS, steady_state

# ---- símbolos ----
P = sp.symbols("SIGG DELTA RSTAR ALFA OMEGA DBAR PSSI PHI STD_EPS_A BETTA RHO")
SIGG, DELTA, RSTAR, ALFA, OMEGA, DBAR, PSSI, PHI, STD_EPS_A, BETTA, RHO = P

_names = ["d", "r", "k", "a", "c", "ivv", "output", "h", "la", "kfu",
          "tb", "tby", "ca", "cay", "tb_o_tot"]
cur = {n: sp.Symbol(n) for n in _names}
nxt = {n: sp.Symbol(n + "p") for n in _names}
g = {**cur, **{n + "p": nxt[n] for n in _names}}   # acesso por nome

# atalhos
d, r, k, a = cur["d"], cur["r"], cur["k"], cur["a"]
c, ivv, output, h, la, kfu = (cur["c"], cur["ivv"], cur["output"], cur["h"],
                              cur["la"], cur["kfu"])
tb, tby, ca, cay, tb_o_tot = (cur["tb"], cur["tby"], cur["ca"], cur["cay"],
                              cur["tb_o_tot"])
dp, rp, kp, ap = nxt["d"], nxt["r"], nxt["k"], nxt["a"]
cp, ivvp, outputp, hp, lap, kfup = (nxt["c"], nxt["ivv"], nxt["output"],
                                    nxt["h"], nxt["la"], nxt["kfu"])
tbp, tbyp, cap, cayp, tb_o_totp = (nxt["tb"], nxt["tby"], nxt["ca"],
                                   nxt["cay"], nxt["tb_o_tot"])

# ---- 15 equações (edeir_model.m) ----
f = sp.Matrix([
    -dp + (1 + r) * d + c + ivv + PHI / 2 * (kp - k) ** 2 - output,          # e1
    -output + a * k ** ALFA * h ** (1 - ALFA),                                # e2
    -la + (c - h ** OMEGA / OMEGA) ** (-SIGG),                                # e3
    -h ** (OMEGA - 1) + (1 - ALFA) * a * (k / h) ** ALFA,                     # e4
    -la + BETTA * (1 + rp) * lap,                                             # e5
    -la * (1 + PHI * (kp - k)) + BETTA * lap * (1 - DELTA
        + ALFA * ap * (kp / hp) ** (ALFA - 1) + PHI * (kfup - kp)),           # e6
    -rp + RSTAR + PSSI * (sp.exp(dp - DBAR) - 1),                             # e7
    -ivv + kp - (1 - DELTA) * k,                                             # e8
    -tb + output - c - ivv,                                                   # e9
    -tby + tb / output,                                                       # e10
    -ca - dp + d,                                                             # e11
    -cay + ca / output,                                                       # e12
    -sp.log(ap) + RHO * sp.log(a),                                            # e13
    -kfu + kp,                                                                # e14
    -tb_o_tot + tb / a,                                                       # e15
])

# ordenação e flags de "em log"
STATE = ["d", "r", "k", "a"]
STATE_LOG = {"d": False, "r": False, "k": True, "a": True}
CONTROL = ["c", "ivv", "output", "h", "la", "kfu",
           "tb", "tby", "ca", "cay", "tb_o_tot"]
CONTROL_LOG = {n: (n in {"c", "ivv", "output", "h", "la", "kfu"}) for n in CONTROL}


def _deriv(vars_dict, names, logflags, isp):
    """Matriz de derivadas de f em relação às variáveis dadas (regra do log)."""
    cols = []
    for n in names:
        v = vars_dict[n + "p"] if isp else vars_dict[n]
        col = f.diff(v)
        if logflags[n]:
            col = col * v            # d/d log(v) = v * d/dv
        cols.append(col)
    return sp.Matrix.hstack(*cols)   # 15 x len(names)


# matrizes simbólicas (uma vez)
_FX = _deriv(g, STATE, STATE_LOG, isp=False)
_FXP = _deriv(g, STATE, STATE_LOG, isp=True)
_FY = _deriv(g, CONTROL, CONTROL_LOG, isp=False)
_FYP = _deriv(g, CONTROL, CONTROL_LOG, isp=True)


def _ss_subs():
    ss = steady_state()
    subs = {}
    for n in _names:
        subs[cur[n]] = ss[n]
        subs[nxt[n]] = ss[n]          # no ss, variável_{t+1} = variável_t
    for pk, pv in PARAMS.items():
        subs[sp.Symbol(pk)] = pv
    return subs


# substitui ss+parâmetros deixando RHO livre; lambdifica em RHO
_SUB = _ss_subs()
_fx_l = sp.lambdify(RHO, _FX.subs(_SUB), "numpy")
_fxp_l = sp.lambdify(RHO, _FXP.subs(_SUB), "numpy")
_fy_l = sp.lambdify(RHO, _FY.subs(_SUB), "numpy")
_fyp_l = sp.lambdify(RHO, _FYP.subs(_SUB), "numpy")
_f_l = sp.lambdify(RHO, f.subs(_SUB), "numpy")


def num_eval(rho):
    """Derivadas numéricas no ss para um dado RHO. Retorna fx, fxp, fy, fyp, eta, fss."""
    fx = np.array(_fx_l(rho), float)
    fxp = np.array(_fxp_l(rho), float)
    fy = np.array(_fy_l(rho), float)
    fyp = np.array(_fyp_l(rho), float)
    eta = np.zeros((4, 1)); eta[3, 0] = PARAMS["STD_EPS_A"]   # choque em a (tot)
    fss = np.array(_f_l(rho), float).ravel()
    return fx, fxp, fy, fyp, eta, fss


if __name__ == "__main__":
    fx, fxp, fy, fyp, eta, fss = num_eval(0.42)
    print("Dimensões: fx", fx.shape, "fxp", fxp.shape, "fy", fy.shape, "fyp", fyp.shape)
    print("Resíduo do estado estacionário max|f(ss)| =", np.max(np.abs(fss)))
    print("fx[0,0] (e1 d/d d) =", fx[0, 0], "(esperado r+1 =", 1 + steady_state()['r'], ")")
    print("fx[12,3] (e13 d/d log a) =", fx[12, 3], "(esperado RHO=0.42)")
    print("fxp[12,3] (e13 d/d log ap) =", fxp[12, 3], "(esperado -1)")
