"""
Camada 0 — Detrending.

Três métodos de remoção de tendência, aplicados de forma uniforme para montar
o dataset da extensão. Cada série-base vira 3 versões (coluna detrending_method).

Regras por variável (seguindo Uribe & Schmitt-Grohé):
  - y, c, i, tot -> tira LOG e depois remove a tendência
  - tb           -> em NÍVEL (pode ser negativo; não se tira log)

Métodos:
  - "quadratic" : remove tendência polinomial de grau 2   (método do livro)
  - "hp"        : filtro Hodrick-Prescott (lambda=100 p/ anual, Ravn-Uhlig)
  - "hamilton"  : filtro de regressão de Hamilton (2018), h=2, p=2 p/ anual

IMPORTANTE sobre a janela: o Uribe detrendou tot em 1980-2011, mas gdp/c/i/tb
numa janela mais longa (1965-2010) e recortou. Portanto, recalcular aqui NÃO
reproduz os números originais dele — para replicação fiel use o usg_tot_data.mat
direto. Este módulo serve à EXTENSÃO (construção nova, consistente).
"""

import numpy as np

# variáveis que entram em log antes do detrend (tb fica em nível)
LOG_VARS = {"tot", "gdp", "c", "i"}
LEVEL_VARS = {"tb"}


# ----------------------------------------------------------------------------
# Métodos (recebem um vetor 1D já transformado -> devolvem o componente cíclico)
# ----------------------------------------------------------------------------
def polynomial_detrend(y, degree=2):
    """Remove tendência polinomial (grau 2 = quadrático, método do livro)."""
    y = np.asarray(y, float)
    n = len(y)
    t = np.arange(n)
    X = np.vander(t, degree + 1, increasing=True)  # [1, t, t^2, ...]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def hp_detrend(y, lamb=100.0):
    """Componente cíclico do filtro HP. lambda=100 é o padrão para dados anuais."""
    y = np.asarray(y, float)
    n = len(y)
    if n < 3:
        return y - y.mean()
    # matriz de segunda diferença D (n-2 x n); minimiza ciclo^2 + lambda*||D*trend||^2
    D = np.zeros((n - 2, n))
    for i in range(n - 2):
        D[i, i], D[i, i + 1], D[i, i + 2] = 1.0, -2.0, 1.0
    trend = np.linalg.solve(np.eye(n) + lamb * D.T @ D, y)
    return y - trend


def hamilton_detrend(y, h=2, p=2):
    """Filtro de regressão de Hamilton (2018).

    Regride y_{t} sobre constante + {y_{t-h}, ..., y_{t-h-p+1}} e devolve o
    resíduo (componente cíclico). Para dados anuais Hamilton sugere h=2, p=2.
    O vetor devolvido tem NaN nas primeiras (h+p-1) posições (sem regressor).
    """
    y = np.asarray(y, float)
    n = len(y)
    start = h + p - 1
    cyc = np.full(n, np.nan)
    if n <= start:
        return cyc
    rows = range(start, n)
    X = np.column_stack([np.ones(len(rows))] +
                        [[y[t - h - j] for t in rows] for j in range(p)])
    yy = np.array([y[t] for t in rows])
    beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
    cyc[start:] = yy - X @ beta
    return cyc


_METHODS = {
    "quadratic": polynomial_detrend,
    "hp": hp_detrend,
    "hamilton": hamilton_detrend,
}


# ----------------------------------------------------------------------------
# Dispatcher: aplica a regra de log/nível + o método escolhido
# ----------------------------------------------------------------------------
def detrend_series(values, varname, method):
    """Detrenda um vetor `values` da variável `varname` pelo `method`.

    Aplica log se varname em LOG_VARS; nível caso contrário.
    """
    if method not in _METHODS:
        raise ValueError(f"método desconhecido: {method!r} (use {list(_METHODS)})")
    x = np.asarray(values, float)
    if varname in LOG_VARS:
        if np.any(x <= 0):
            raise ValueError(f"{varname}: valores <=0 impedem log")
        x = np.log(x)
    return _METHODS[method](x)


if __name__ == "__main__":
    # teste de sanidade: série log-quadrática pura -> quadratic recupera o ciclo
    rng = np.random.default_rng(0)
    n = 32
    t = np.arange(n)
    cycle = rng.normal(size=n) * 0.1
    log_level = 5 + 0.03 * t - 0.0005 * t**2 + cycle   # tendência quadrática + ciclo
    level = np.exp(log_level)

    q = detrend_series(level, "gdp", "quadratic")
    print("quadratic recupera ciclo? corr =",
          round(float(np.corrcoef(q, cycle)[0, 1]), 4),
          "| max erro =", round(float(np.max(np.abs(q - cycle))), 4))

    hp = detrend_series(level, "gdp", "hp")
    ham = detrend_series(level, "gdp", "hamilton")
    print("HP        -> desvio padrão do ciclo:", round(float(np.std(hp)), 4))
    print("Hamilton  -> n válidos:", int(np.sum(~np.isnan(ham))), "de", n)

    # tb em nível (pode ser negativo)
    tb = np.array([-0.02, 0.01, 0.03, -0.01, 0.0, 0.02, -0.03, 0.01] * 4)[:n]
    tbq = detrend_series(tb, "tb", "quadratic")
    print("tb (nível) quadratic ok? média ~0:", round(float(tbq.mean()), 6))
