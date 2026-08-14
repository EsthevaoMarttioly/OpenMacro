"""
Camada dsge — Solver de 1ª ordem gx_hx (perturbação / método de Klein 2000).

Tradução de gx_hx.m (Schmitt-Grohé & Uribe), substituindo a rotina QZ do MATLAB
por scipy.linalg.ordqz. Resolve o sistema log-linearizado
    E_t f(yp, y, xp, x) = 0
para a solução na forma estado-espaço:
    y_t     = gx x_t
    x_{t+1} = hx x_t  (+ choques)
onde x são os estados e y os controles.

Baseia-se na condição  fx + fxp*hx + fy*gx + fyp*gx*hx = 0, resolvida via
decomposição de Schur generalizada ordenada (autovalores estáveis primeiro).
"""

import numpy as np
from scipy.linalg import ordqz


def gx_hx(fy, fx, fyp, fxp, stake=1.0):
    fy = np.asarray(fy, float); fx = np.asarray(fx, float)
    fyp = np.asarray(fyp, float); fxp = np.asarray(fxp, float)
    nk = fx.shape[1]                      # número de estados

    # Sistema A0 w_{t+1} = B0 w_t, com w = [x; y]
    A0 = np.hstack([fxp, fyp])
    B0 = -np.hstack([fx, fy])

    # A dinâmica é w_{t+1} = A0^{-1} B0 w_t. Seus autovalores mu satisfazem
    # B0 z = mu A0 z. Em termos do par (A0, B0), lambda = alpha/beta e mu = 1/lambda,
    # logo |mu| < stake  <=>  |beta| < stake*|alpha|. Ordenamos estáveis primeiro.
    def sort(alpha, beta):
        return np.abs(beta) < stake * np.abs(alpha)

    S, T, alpha, beta, Q, Z = ordqz(A0, B0, sort=sort, output="complex")

    n_stable = int(np.sum(np.abs(beta) < stake * np.abs(alpha)))
    exitflag = 1 if n_stable == nk else 0   # 1 = existência+unicidade (Blanchard-Kahn)

    z11 = Z[:nk, :nk]
    z21 = Z[nk:, :nk]
    s11 = S[:nk, :nk]
    t11 = T[:nk, :nk]

    z11i = np.linalg.inv(z11)
    gx = np.real(z21 @ z11i)
    hx = np.real(z11 @ np.linalg.solve(s11, t11) @ z11i)
    return gx, hx, exitflag


def residual(gx, hx, fy, fx, fyp, fxp):
    """Resíduo da condição de equilíbrio: deve ser ~0 se a solução está correta."""
    return fx + fxp @ hx + fy @ gx + fyp @ gx @ hx


if __name__ == "__main__":
    from edeir_model import num_eval
    fx, fxp, fy, fyp, eta, _ = num_eval(0.42)
    gx, hx, flag = gx_hx(fy, fx, fyp, fxp)
    R = residual(gx, hx, fy, fx, fyp, fxp)
    print("exitflag (1=ok, Blanchard-Kahn):", flag)
    print("nº estados:", hx.shape[0], "| resíduo max|f-cond| =", np.max(np.abs(R)))
    print("autovalores de hx (devem ter |.|<1):", np.round(np.abs(np.linalg.eigvals(hx)), 4))
