"""
mom — momentos de segunda ordem incondicionais.

Tradução de mom.m (Schmitt-Grohé & Uribe). Para o sistema
    x_{t+1} = hx x_t + e_{t+1},   E[e e'] = varshock
    y_t     = gx x_t
computa a matriz de covariância incondicional de x(t) com x(t+J), isto é
sigxJ = E[x_t x_{t+J}'], e a de y, sigyJ = E[y_t y_{t+J}'].

Usa o "doubling algorithm" (padrão no MATLAB do Uribe), que resolve a equação de
Lyapunov  sigx = hx sigx hx' + varshock  de forma rápida e estável.

Função pura: não lê dados.
"""

import numpy as np


def mom(gx, hx, varshock, J=0, method=1, tol=1e-25, maxiter=1000):
    gx = np.atleast_2d(np.asarray(gx, float))
    hx = np.atleast_2d(np.asarray(hx, float))
    varshock = np.atleast_2d(np.asarray(varshock, float))

    if method == 1:
        # ----- doubling algorithm -----
        hx_old = hx.copy()
        sig_old = varshock.copy()
        sigx_old = np.eye(hx.shape[0])
        diff = 0.1
        it = 0
        while diff > tol and it < maxiter:
            sigx = hx_old @ sigx_old @ hx_old.T + sig_old
            diff = np.max(np.abs(sigx - sigx_old))
            sig_old = hx_old @ sig_old @ hx_old.T + sig_old
            hx_old = hx_old @ hx_old
            sigx_old = sigx
            it += 1
    else:
        # ----- método algébrico (Kronecker) -----
        n = hx.shape[0]
        F = np.kron(hx, hx)
        vec = np.linalg.solve(np.eye(F.shape[0]) - F, varshock.reshape(-1, order="F"))
        sigx = vec.reshape(n, n, order="F")

    if J != 0:
        # E[x_t x_{t+J}'] = hx^(-min(0,J)) * sigx * (hx')^max(0,J)
        left = np.linalg.matrix_power(hx, -min(0, J))
        right = np.linalg.matrix_power(hx.T, max(0, J))
        sigxJ = left @ sigx @ right
    else:
        sigxJ = sigx

    sigyJ = np.real(gx @ sigxJ @ gx.T)
    sigxJ = np.real(sigxJ)
    return sigyJ, sigxJ


if __name__ == "__main__":
    # AR(1) escalar: var incondicional = sigma^2 / (1 - rho^2)
    rho, s2 = 0.7, 1.0
    _, sigx = mom(np.array([[1.0]]), np.array([[rho]]), np.array([[s2]]))
    print("sigx:", sigx.item(), "| esperado:", s2 / (1 - rho**2))
    assert abs(sigx.item() - s2 / (1 - rho**2)) < 1e-10

    # confere doubling == kronecker num caso 2x2
    hx = np.array([[0.5, 0.1], [0.0, 0.3]])
    vs = np.array([[1.0, 0.2], [0.2, 0.5]])
    _, s_dbl = mom(np.eye(2), hx, vs, method=1)
    _, s_krn = mom(np.eye(2), hx, vs, method=2)
    assert np.allclose(s_dbl, s_krn, atol=1e-10)
    print("OK")
