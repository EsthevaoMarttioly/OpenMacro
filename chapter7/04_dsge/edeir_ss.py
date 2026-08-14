"""
Camada dsge — Calibração e estado estacionário do modelo SOE-RBC (EDEIR).

Tradução de edeir_ss.m (Uribe). Modelo de economia aberta com juro elástico à
dívida (External Debt-Elastic Interest Rate) do cap. 4, com o choque de
produtividade `a` reinterpretado como termos de troca (tot). Calibração anual.
"""

import numpy as np

# ---- Calibração (edeir_ss.m) ----
SIGG = 2.0            # aversão ao risco (Mendoza)
DELTA = 0.1          # depreciação
RSTAR = 0.04         # juro de longo prazo
ALFA = 0.32          # F(k,h) = k^ALFA h^(1-ALFA)
OMEGA = 1.455        # Frisch (Mendoza 1991)
DBAR = 0.74421765717098   # dívida
PSSI = 0.11135 / 150      # elasticidade do juro à dívida
PHI = 0.028          # custo de ajuste do capital
STD_EPS_A = 0.0129   # desvio-padrão da inovação do choque
BETTA = 1.0 / (1.0 + RSTAR)

PARAMS = dict(SIGG=SIGG, DELTA=DELTA, RSTAR=RSTAR, ALFA=ALFA, OMEGA=OMEGA,
              DBAR=DBAR, PSSI=PSSI, PHI=PHI, STD_EPS_A=STD_EPS_A, BETTA=BETTA)


def steady_state():
    r = RSTAR
    d = DBAR
    KAPA = ((1.0 / BETTA - (1.0 - DELTA)) / ALFA) ** (1.0 / (ALFA - 1.0))
    h = ((1.0 - ALFA) * KAPA ** ALFA) ** (1.0 / (OMEGA - 1.0))
    k = KAPA * h
    output = KAPA ** ALFA * h
    c = output - DELTA * k - RSTAR * DBAR
    ivv = DELTA * k
    tb = output - ivv - c
    tby = tb / output
    tb_o_tot = tb
    ca = -r * d + tb
    cay = ca / output
    a = 1.0
    la = (c - h ** OMEGA / OMEGA) ** (-SIGG)
    kfu = k
    return dict(d=d, r=r, k=k, a=a, c=c, ivv=ivv, output=output, h=h, la=la,
                kfu=kfu, tb=tb, tby=tby, ca=ca, cay=cay, tb_o_tot=tb_o_tot)


if __name__ == "__main__":
    ss = steady_state()
    for kk, vv in ss.items():
        print(f"  {kk:9s} = {vv: .6f}")
