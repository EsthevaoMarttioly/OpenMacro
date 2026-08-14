"""
Camada 05_mx — Tabelas 7.6 e 7.7: variância explicada por choques de tot,
modelo MX vs SVAR empírico.

Tradução de tot_mx_var_cbc.m (Uribe). Para cada país resolve o MX model com a
calibração de phi_pssi_cbc.mat e compara a variância condicional a tot (do modelo)
com a variância total e a fatia de tot do SVAR empírico (tot_svar_var_cbc.mat).

Linhas: [tot, tby, output, c, ivv]  (= [tot, tb, y, c, i]).
"""

from pathlib import Path
import numpy as np
import scipy.io as sio

from mx_function import moments
from mx_model import IDX_STATE, IDX_CONTROL

MX = Path(__file__).resolve().parents[2] / "usg_mx"
PKG = Path(__file__).resolve().parents[1]

# variáveis do modelo (constant prices), na ordem [tot, tby, output, c, ivv]
ROWS = ["tby_constant_prices", "output_constant_prices",
        "c_constant_prices", "ivv_constant_prices"]
LABELS = ["tot", "tb", "y", "c", "i"]


def build():
    phi = sio.loadmat(MX / "phi_pssi_cbc.mat")
    PHIM = phi["PHIM_cbc"].ravel(); PHIX = phi["PHIX_cbc"].ravel()
    PSSI = phi["PSSI_cbc"].ravel(); RHO = phi["RHO_TOT_cbc"].ravel()
    STD = phi["STD_TOT_cbc"].ravel()

    sv = sio.loadmat(MX / "tot_svar_var_cbc.mat")
    var_svar = sv["var_svar"]                      # (5, ncou) variância total empírica
    svar_ratio = sv["var_svar_tot_ratio_cbc"]      # (5, ncou) fatia de tot no SVAR
    ncou = int(sv["ncou"].ravel()[0])
    country = [str(sv["country"][0, i][0]) for i in range(sv["country"].shape[1])]
    i_c = sv["i_c"].ravel().astype(int) - 1        # índice (base 0)

    ntot = IDX_STATE["tot"]
    var_model = np.zeros((5, ncou))
    for k in range(ncou):
        sigy0, sigx0, flag = moments(PHIM[k], PHIX[k], PSSI[k], RHO[k], STD[k])
        var_model[0, k] = sigx0[ntot, ntot]
        for j, nm in enumerate(ROWS, start=1):
            var_model[j, k] = sigy0[IDX_CONTROL[nm], IDX_CONTROL[nm]]

    model_ratio = var_model / var_svar             # fração da var. empírica total
    tbl76_model = np.median(model_ratio, axis=1) * 100
    tbl76_svar = np.median(svar_ratio, axis=1) * 100
    names = [country[i_c[k]] for k in range(ncou)]
    return tbl76_model, tbl76_svar, model_ratio, svar_ratio, names


if __name__ == "__main__":
    m76, s76, mr, sr, names = build()
    print("=== Tabela 7.6 (medianas, %) ===")
    print(f"{'':6s}{'MX model':>12s}{'SVAR':>10s}")
    for i, l in enumerate(LABELS):
        print(f"{l:6s}{m76[i]:12.1f}{s76[i]:10.1f}")

    outdir = PKG / "05_mx" / "results"
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "table_7_6.txt", "w") as f:
        f.write(f"{'':6s}{'MX model':>12s}{'SVAR':>10s}\n")
        for i, l in enumerate(LABELS):
            f.write(f"{l:6s}{m76[i]:12.1f}{s76[i]:10.1f}\n")

    import csv
    with open(outdir / "table_7_7_bycountry.csv", "w", newline="") as f:
        wr = csv.writer(f)
        hdr = ["country"]
        for l in LABELS:
            hdr += [f"{l}_model", f"{l}_svar"]
        wr.writerow(hdr)
        for k, nm in enumerate(names):
            row = [nm]
            for i in range(5):
                row += [f"{mr[i,k]*100:.1f}", f"{sr[i,k]*100:.1f}"]
            wr.writerow(row)
    print(f"\nSalvo em {outdir.relative_to(PKG)}")
