"""
lagg — monta série contemporânea + defasagens.

Tradução de lagg.m (Uribe). Dado y de dimensão T x v (primeira linha = data mais
antiga) e L >= 1 inteiro, devolve uma matriz (T-L) x (v*(L+1)) onde as primeiras
v colunas são a série contemporânea, as v seguintes são a série defasada 1
período, e assim por diante.

É uma função pura: não lê dados, só transforma o array que recebe.
"""

import numpy as np


def lagg(y, L=1):
    y = np.asarray(y, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    T, v = y.shape
    # j=0 -> contemporâneo (linhas L..T-1); j=k -> defasagem k
    blocks = [y[L - j: T - j, :] for j in range(L + 1)]
    return np.hstack(blocks)


if __name__ == "__main__":
    y = np.arange(1, 7).reshape(-1, 1).astype(float)  # [1;2;3;4;5;6]
    out = lagg(y, 1)
    print("y:", y.ravel())
    print("lagg(y,1):\n", out)
    # esperado: col0 = [2,3,4,5,6] (contemp), col1 = [1,2,3,4,5] (lag1)
    assert np.allclose(out[:, 0], [2, 3, 4, 5, 6])
    assert np.allclose(out[:, 1], [1, 2, 3, 4, 5])
    # L=2
    out2 = lagg(y, 2)
    assert out2.shape == (4, 3)
    assert np.allclose(out2[:, 0], [3, 4, 5, 6])   # contemp
    assert np.allclose(out2[:, 1], [2, 3, 4, 5])   # lag1
    assert np.allclose(out2[:, 2], [1, 2, 3, 4])   # lag2
    print("OK")
