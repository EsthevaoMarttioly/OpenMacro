"""
ir — funções de resposta a impulso (impulse responses).

Tradução de ir.m (Schmitt-Grohé & Uribe). Computa as respostas de T períodos do
vetor [y; x], cuja lei de movimento é:
    x_{t+1} = hx x_t
    y_t     = gx x_t
com condição inicial x0.

Devolve:
    IR  : T x (ny+nx), colunas = [y ; x] a cada período (linha t = período t)
    IRy : T x ny  (parte de y)
    IRx : T x nx  (parte de x)

Função pura: não lê dados.
"""

import numpy as np


def ir(gx, hx, x0, T=10):
    gx = np.atleast_2d(np.asarray(gx, float))
    hx = np.atleast_2d(np.asarray(hx, float))
    x0 = np.asarray(x0, float).reshape(-1)
    n = x0.size

    MX = np.vstack([gx, np.eye(n)])          # empilha y (via gx) e x (via I)
    IR = np.zeros((T, MX.shape[0]))
    x = x0.copy()
    for t in range(T):
        IR[t, :] = MX @ x
        x = hx @ x

    IRx = IR[:, -n:]
    IRy = IR[:, :-n]
    return IR, IRy, IRx


if __name__ == "__main__":
    # AR(1) escalar: x_{t+1}=rho x_t, y_t = x_t. Resposta a x0=1 deve ser rho^t.
    rho = 0.5
    IR, IRy, IRx = ir(np.array([[1.0]]), np.array([[rho]]), [1.0], T=5)
    print("IRx (deve ser rho^t):", IRx.ravel())
    assert np.allclose(IRx.ravel(), rho ** np.arange(5))
    print("OK")
