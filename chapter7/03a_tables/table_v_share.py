"""
Camada 3a — Tabela 7.2: participação da variância explicada por choques de tot.

Tradução de table_v_share_cbc.m (Uribe). A partir do SVAR de 5 variáveis
(tot_cbcs sizevar=5) e do painel (totss), reporta a fatia da variância de
[tot, tb, gdp, c, i] explicada por choques de tot, de 3 formas:
  - Cross-Country Mean / Median / MAD das fatias país-a-país
  - Usando a MÉDIA entre países de hx e PI (um SVAR "representativo")
  - Estimação em painel (totss)

quadratic -> 03a_tables/results/ ; hp, hamilton -> extension/tables/
"""

from pathlib import Path
import sys
import numpy as np

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG / "01_primitives"))
from variance_decomposition import variance_decomposition   # noqa: E402

LABELS = ["tot", "tb", "gdp", "c", "i"]


def _dir(method, sub):
    return PKG / (f"02_svar/{sub}" if method == "quadratic" else f"extension/{sub}")


def table_v_share(method="quadratic"):
    # tot_cbcs5
    svar_dir = PKG / ("02_svar/results" if method == "quadratic" else "extension/svar")
    z5 = np.load(svar_dir / f"tot_cbcs_{method}_5.npz", allow_pickle=True)
    vs = z5["v_share"]                       # (ncou, 5)
    hx = z5["hx"]                            # (ncou, 5, 5)
    PI = z5["PI"]

    mean_ = vs.mean(axis=0)
    median_ = np.median(vs, axis=0)
    # O livro rotula "Median Absolute Deviation", mas o número vem do mad() do
    # MATLAB, que por padrão é o desvio absoluto MÉDIO (sobre a média).
    # Usamos essa definição para reproduzir a Tabela 7.2 exatamente.
    mad_ = np.mean(np.abs(vs - mean_), axis=0)

    # SVAR "representativo": média entre países de hx e PI
    H = hx.mean(axis=0)
    P = PI.mean(axis=0)
    _, vdx, _, _ = variance_decomposition(np.eye(5), H, P)
    mean_hxpi = vdx[0, :] * 100.0

    # painel
    zt = np.load(svar_dir / f"totss_{method}.npz", allow_pickle=True)
    panel = zt["v_share"]

    return {
        "Cross-Country Mean": mean_,
        "Cross-Country Median": median_,
        "Median Absolute Deviation": mad_,
        "Using Mean of hx and PI": mean_hxpi,
        "Panel Estimation": panel,
    }, vs, z5["countries"]


def _fmt(rows):
    lines = [f"{'':28s}" + "".join(f"{l:>7s}" for l in LABELS)]
    for name, vals in rows.items():
        lines.append(f"{name:28s}" + "".join(f"{v:7.1f}" for v in vals))
    return "\n".join(lines)


if __name__ == "__main__":
    import csv
    for method in ["quadratic", "hp", "hamilton"]:
        rows, vs, countries = table_v_share(method)
        outdir = PKG / ("03a_tables/results" if method == "quadratic"
                        else "extension/tables")
        outdir.mkdir(parents=True, exist_ok=True)
        txt = _fmt(rows)
        (outdir / f"table_7_2_vshare_{method}.txt").write_text(txt + "\n")
        with open(outdir / f"table_7_2_vshare_{method}_bycountry.csv", "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["country"] + LABELS)
            for c, row in zip(countries, vs):
                wr.writerow([c] + [f"{x:.4f}" for x in row])
        print(f"=== Tabela 7.2  [{method}] ===")
        print(txt, "\n")
