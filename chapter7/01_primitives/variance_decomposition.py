"""
variance_decomposition — decomposição de variância.

Tradução de variance_decomposition.m (Uribe). Para o sistema
    x_{t+1} = hx x_t + ETA1 eps_{t+1},   eps ~ iid N(0, I)
    y_t     = gx x_t                     (ou y_t = gx x_t + ETA2 mu_t, opcional)
computa quanto cada choque (coluna de ETA1) explica da variância de cada
variável de x e de y.

Retorna:
    Vyr, Vxr : decomposições em FRAÇÃO da variância total
               elemento (i, j) = parcela da variância da variável j explicada
               pelo choque i.
    Vy,  Vx  : decomposições em NÍVEL (variâncias)

Função pura: não lê dados (usa mom internamente).
"""

import numpy as np

try:
    from .mom import mom
except ImportError:
    from mom import mom


def variance_decomposition(gx, hx, ETA1, ETA2=None):
    gx = np.atleast_2d(np.asarray(gx, float))
    hx = np.atleast_2d(np.asarray(hx, float))
    ETA1 = np.atleast_2d(np.asarray(ETA1, float))

    n1 = ETA1.shape[1]
    Vy, Vx = [], []

    for j in range(n1):
        I1 = np.zeros((n1, n1))
        I1[j, j] = 1.0
        V1 = ETA1 @ I1 @ ETA1.T            # covariância só do choque j
        sigy, sigx = mom(gx, hx, V1)
        Vy.append(np.diag(sigy).copy())
        Vx.append(np.diag(sigx).copy())

    Vy = np.array(Vy)
    Vx = np.array(Vx)

    if ETA2 is not None:
        ETA2 = np.atleast_2d(np.asarray(ETA2, float))
        n2 = ETA2.shape[1]
        extra = []
        for j in range(n2):
            I2 = np.zeros((n2, n2))
            I2[j, j] = 1.0
            V2 = ETA2 @ I2 @ ETA2.T
            extra.append(np.diag(V2).copy())
        Vy = np.vstack([Vy, np.array(extra)])

    Vyr = Vy / Vy.sum(axis=0, keepdims=True)
    Vxr = Vx / Vx.sum(axis=0, keepdims=True)
    return Vyr, Vxr, Vy, Vx


if __name__ == "__main__":
    # dois estados independentes, um choque cada -> cada choque explica 100%
    # da sua própria variável e 0% da outra.
    hx = np.array([[0.5, 0.0], [0.0, 0.3]])
    ETA1 = np.array([[1.0, 0.0], [0.0, 1.0]])
    Vyr, Vxr, Vy, Vx = variance_decomposition(np.eye(2), hx, ETA1)
    print("Vxr (shares):\n", np.round(Vxr, 6))
    assert np.allclose(Vxr, np.eye(2), atol=1e-8)
    assert np.allclose(Vxr.sum(axis=0), 1.0)
    print("OK")
