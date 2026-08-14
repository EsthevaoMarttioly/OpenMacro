"""
Camada 05_mx — Estado estacionário do MX Model (dois setores).

Tradução de mx_ss.m (Uribe). Modelo de economia aberta com um bem importável (m)
e um exportável (x), agregador Armington Cobb-Douglas (MU=1), juro elástico à
dívida, dirigido por choques de termos de troca (tot).

Os parâmetros PHIM, PHIX, PSSI, RHO_TOT, STD_TOT são calibrados por país
(phi_pssi_cbc.mat); os demais são fixos (calibração do livro).
"""

import numpy as np

# ---- parâmetros fixos (mx_ss.m) ----
SIGG = 2.0
OMEGA = 1.455
ALFA = 0.32
MU = 1.0            # Cobb-Douglas (Armington)
RSTAR = 0.11
TFP_M = 1.0
TOTBAR = 1.0
DELTA = 0.1
SX = 0.26 * 0.81    # participação das exportações no PIB
STB = 0.01          # razão balança comercial / produto no ss
GAMA = 0.88         # PIB do setor X sobre PIB do setor M
BETTA = 1.0 / (1.0 + RSTAR)

# medianas cross-country (usadas quando não se passa calibração de país)
PHIM_MED, PHIX_MED, PSSI_MED, RHO_TOT_MED, STD_TOT_MED = 1.82, 1.56, 0.18, 0.5, 0.1


def steady_state(PHIM=PHIM_MED, PHIX=PHIX_MED, PSSI=PSSI_MED,
                 RHO_TOT=RHO_TOT_MED, STD_TOT=STD_TOT_MED):
    tot = TOTBAR
    um = 1.0 / BETTA - 1.0 + DELTA
    ux = 1.0 / BETTA - 1.0 + DELTA

    A = GAMA ** ((1 - ALFA) * (OMEGA - 1) / OMEGA) / TOTBAR
    tfp_m = TFP_M
    tfp_x = A * tfp_m
    TFP_X = tfp_x

    UPSILON = (GAMA / (1 + GAMA) - SX) / (1.0 / (1 + GAMA) + SX - STB)
    UPSILON_tot = (UPSILON * tot ** (MU - 1)) ** (-1.0 / MU)
    CHI = UPSILON_tot / (1 + UPSILON_tot)

    ax_o_am = (CHI / (1 - CHI) * tot) ** (-MU)
    if MU == 1:
        pm = CHI * ax_o_am ** (1 - CHI)
    else:
        pm = (CHI + (1 - CHI) * ax_o_am ** (1 - 1 / MU)) ** (1 / (1 - 1 / MU) - 1) * CHI
    px = tot * pm

    km_o_hm = (um / pm / ALFA / tfp_m) ** (1 / (ALFA - 1))
    wm = pm * (1 - ALFA) * tfp_m * km_o_hm ** ALFA
    kx_o_hx = (ux / px / ALFA / tfp_x) ** (1 / (ALFA - 1))
    wx = tfp_x * px * (1 - ALFA) * kx_o_hx ** ALFA

    hm = wm ** (1 / (OMEGA - 1))
    hx = wx ** (1 / (OMEGA - 1))
    kx = hx * kx_o_hx
    km = hm * km_o_hm
    ym = tfp_m * km ** ALFA * hm ** (1 - ALFA)
    yx = tfp_x * kx ** ALFA * hx ** (1 - ALFA)
    im = DELTA * km
    ix = DELTA * kx
    ivv = im + ix
    y = px * yx + pm * ym
    r = RSTAR
    d = STB * y * (1 + r) / r
    DBAR = d
    m = y * (SX - STB)
    am = ym + m / pm
    x = SX * y
    ax = yx - x / px
    if MU == 1:
        c = am ** CHI * ax ** (1 - CHI) - im - ix
    else:
        c = (CHI * am ** (1 - 1 / MU) + (1 - CHI) * ax ** (1 - 1 / MU)) ** (1 / (1 - 1 / MU)) - im - ix
    la = (c - hm ** OMEGA / OMEGA - hx ** OMEGA / OMEGA) ** (-SIGG)
    output = y
    tby = (x - m) / output

    ss = dict(
        d=d, km=km, kx=kx, tfp_m=tfp_m, tfp_x=tfp_x, tot=tot,
        c=c, hm=hm, h_x=hx, im=im, ix=ix, am=am, ax=ax, px=px, pm=pm,
        ym=ym, yx=yx, wm=wm, wx=wx, ux=ux, um=um, m=m, x=x, output=output,
        la=la, kmfu=km, kxfu=kx, ivv=ivv,
        output_constant_prices=output, c_constant_prices=c,
        ivv_constant_prices=ivv, tby=tby, r=r, tby_constant_prices=tby,
    )
    consts = dict(SIGG=SIGG, OMEGAM=OMEGA, OMEGAX=OMEGA, ALFAM=ALFA, ALFAX=ALFA,
                  BETTA=BETTA, DELTA=DELTA, CHI=CHI, MU=MU, RSTAR=RSTAR,
                  DBAR=DBAR, PSSI=PSSI, PHIM=PHIM, PHIX=PHIX,
                  RHO_TOT=RHO_TOT, STD_TOT=STD_TOT,
                  TOT=TOTBAR, TFP_X=TFP_X, TFP_M=TFP_M, OUTPUT=output,
                  PX=px, PM=pm)
    return ss, consts


if __name__ == "__main__":
    ss, consts = steady_state()
    print("CHI =", round(consts["CHI"], 5), "| tby =", round(ss["tby"], 5),
          "| output =", round(ss["output"], 5))
    for k in ["c", "km", "kx", "hm", "h_x", "am", "ax", "px", "pm", "la", "d", "m", "x"]:
        print(f"  {k:8s} = {ss[k]: .6f}")
