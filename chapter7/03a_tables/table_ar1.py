"""
Camada 3a — Tabela 7.1: processo AR(1) dos termos de troca.

Tradução de table_ar1_cbc.m (Uribe). A partir das estimativas país-a-país do
SVAR de 2 variáveis (tot_cbcs sizevar=2), resume o processo do tot:
    rho = hx[0,0]   (persistência)
    pi  = PI[0,0]   (desvio-padrão da inovação)
    R2  = R2[0]     (ajuste)
reportando média, mediana e intervalo interquartílico (mesmo método do MATLAB:
valor nas posições round(ncou*0.25) e round(ncou*0.75) da série ordenada).

quadratic -> 03a_tables/results/ ; hp, hamilton -> extension/tables/
"""

from pathlib import Path
import numpy as np

PKG = Path(__file__).resolve().parents[1]


def _npz(method, sizevar):
    d = PKG / ("02_svar/results" if method == "quadratic" else "extension/svar")
    return np.load(d / f"tot_cbcs_{method}_{sizevar}.npz", allow_pickle=True)


def _iqr(x):
    """IQR no estilo do MATLAB: valores nas posições round(n*0.25) e round(n*0.75)."""
    xs = np.sort(x)
    n = len(x)
    q25 = xs[int(np.round(n * 0.25)) - 1]
    q75 = xs[int(np.round(n * 0.75)) - 1]
    return q25, q75


def table_ar1(method="quadratic"):
    z = _npz(method, 2)
    rho = z["hx"][:, 0, 0]
    pi = z["PI"][:, 0, 0]
    r2 = z["R2"][:, 0]

    rows = {}
    for name, x in [("rho", rho), ("pi", pi), ("R2", r2)]:
        q25, q75 = _iqr(x)
        rows[name] = dict(mean=x.mean(), median=np.median(x), q25=q25, q75=q75)
    return rows, dict(rho=rho, pi=pi, R2=r2, countries=z["countries"])


def _fmt(rows):
    lines = []
    lines.append(f"{'':22s}{'rho':>10s}{'pi':>10s}{'R2':>10s}")
    lines.append(f"{'Mean':22s}" + "".join(f"{rows[k]['mean']:10.2f}" for k in ['rho','pi','R2']))
    lines.append(f"{'Median':22s}" + "".join(f"{rows[k]['median']:10.2f}" for k in ['rho','pi','R2']))
    lines.append(f"{'IQR (25th pct)':22s}" + "".join(f"{rows[k]['q25']:10.2f}" for k in ['rho','pi','R2']))
    lines.append(f"{'IQR (75th pct)':22s}" + "".join(f"{rows[k]['q75']:10.2f}" for k in ['rho','pi','R2']))
    return "\n".join(lines)


if __name__ == "__main__":
    for method in ["quadratic", "hp", "hamilton"]:
        rows, detail = table_ar1(method)
        outdir = PKG / ("03a_tables/results" if method == "quadratic"
                        else "extension/tables")
        outdir.mkdir(parents=True, exist_ok=True)
        txt = _fmt(rows)
        (outdir / f"table_7_1_ar1_{method}.txt").write_text(txt + "\n")
        # csv país-a-país
        import csv
        with open(outdir / f"table_7_1_ar1_{method}_bycountry.csv", "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["country", "rho", "pi", "R2"])
            for c, a, b, r in zip(detail["countries"], detail["rho"],
                                  detail["pi"], detail["R2"]):
                wr.writerow([c, f"{a:.6f}", f"{b:.6f}", f"{r:.6f}"])
        print(f"=== Tabela 7.1  [{method}] ===")
        print(txt, "\n")
