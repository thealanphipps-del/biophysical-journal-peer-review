#!/usr/bin/env python3
r"""
Generates the Publication-Ready 4-Panel Headline Figure (Figure 1) for the N = 1,040 Cohort
Biophysical Journal (Cell Press) White-Background Vector & High-Resolution Print Standard
(Release ID: SOV-BIOPHYS-N1040-RC2-20261019)

Compliance Guarantees:
1. Pure White Background (#FFFFFF) with Black Axes/Text (#000000), Deep Cyan (#008080) for
   Catalytic Enzymes (N = 922), and Crimson (#D62828) for Structural/Fibrillar Proteins (N = 118).
2. TrueType Font Embedding (pdf.fonttype = 42, ps.fonttype = 42) so all vector/text labels are
   100% searchable and immune to OCR or rasterization artifacts.
3. Explicit Distance-Weighted Kirchhoff Contact Laplacian (dw-GNM) terminology and dimensionless
   scaling prefactor formulation: \nu_1 (THz) = a_0 (\langle z \rangle / (R_g / R_0))^b with R_0 = 1.0 \AA.
4. Exports 600 DPI PNG, Vector PDF, and 600 DPI TIFF.
"""

import json
import shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LEGACY_DIR = Path(r"C:\Users\theal\.gemini\antigravity\brain\337b4327-e578-4669-b7fb-f2396755ea1c")
ACTIVE_DIR = Path(r"C:\Users\theal\.gemini\antigravity\brain\b67c9874-0a7d-48c7-9eb2-fe98a3345a81")

JSON_1000 = LEGACY_DIR / "pdb_1000_cohort_validation.json"
if not JSON_1000.exists():
    JSON_1000 = ACTIVE_DIR / "pdb_1000_cohort_validation.json"


def main():
    data = json.loads(JSON_1000.read_text(encoding="utf-8"))
    records = data["cohort_records"]

    z_rg = np.array([r["compactness_z_over_rg"] for r in records], dtype=np.float64)
    nu_calc = np.array([r["nu_1_thz"] for r in records], dtype=np.float64)
    nu_pred = np.array([r["nu_1_pred_thz"] for r in records], dtype=np.float64)
    is_enz = np.array([r["functional_group"] == "Enzyme" for r in records], dtype=bool)

    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "font.size": 11.5,
        "axes.titlesize": 12.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 11.5,
        "axes.labelweight": "bold",
        "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5,
        "figure.facecolor": "#FFFFFF",
        "axes.facecolor": "#FFFFFF",
        "savefig.facecolor": "#FFFFFF",
        "axes.edgecolor": "#000000",
        "axes.linewidth": 1.1,
        "axes.labelcolor": "#000000",
        "xtick.color": "#000000",
        "ytick.color": "#000000",
        "text.color": "#000000",
        "grid.color": "#CBD5E1",
        "grid.linestyle": "--",
        "grid.alpha": 0.75,
    })

    fig, axes = plt.subplots(2, 2, figsize=(16.0, 12.0), dpi=600)
    fig.suptitle(
        "Three-Dimensional Contact Topology Predicts dw-GNM Kirchhoff Fundamental Frequencies (N = 1,040)\n"
        r"Dimensionless Scaling Law: $\nu_1\ (\mathrm{THz}) = a_0 \left(\langle z \rangle / (R_g / R_0)\right)^{b}$   "
        r"($a_0 = 0.42338\ \mathrm{THz},\ R_0 = 1.0\ \mathrm{\AA},\ b = 1.1194,\ R^2 = 0.8048,\ R^2_{\mathrm{CV}} = 0.8039$)",
        fontsize=13.5,
        fontweight="bold",
        color="#000000",
        y=0.978,
    )

    # =========================================================================
    # Panel A: dw-GNM Kirchhoff Fundamental Frequency vs. Predicted Frequency
    # =========================================================================
    ax1 = axes[0, 0]
    ax1.grid(True)
    ax1.scatter(
        nu_pred[is_enz], nu_calc[is_enz],
        s=26, alpha=0.72, color="#008080", edgecolors="none",
        label=f"Catalytic Enzymes (N = {int(np.sum(is_enz))}, b = 1.16)"
    )
    ax1.scatter(
        nu_pred[~is_enz], nu_calc[~is_enz],
        s=38, alpha=0.90, color="#D62828", marker="^", edgecolors="none",
        label=f"Structural / Fibrillar / Regulatory (N = {int(np.sum(~is_enz))}, b = 0.90)"
    )
    lim_max = max(float(np.max(nu_calc)), float(np.max(nu_pred))) * 1.05
    line_x = np.linspace(0.01, lim_max, 200)
    ax1.plot(line_x, line_x, color="#000000", lw=2.0, ls="--", label="1:1 Identity Line (Calculated = Predicted)")
    ax1.fill_between(
        line_x,
        np.maximum(0, line_x - 0.0354),
        line_x + 0.0354,
        color="#64748B",
        alpha=0.18,
        label=r"$\pm 1\ \mathrm{RMSE\ Band}\ (\pm 0.035\ \mathrm{THz})$"
    )
    ax1.set_title(r"A. dw-GNM Kirchhoff $\nu_1$ vs. Scaling-Law Prediction ($N = 1,040$)", pad=12)
    ax1.set_xlabel(r"Predicted $\nu_1$ from $a_0 \left(\langle z \rangle / (R_g / R_0)\right)^{1.12}$ (THz)", labelpad=7)
    ax1.set_ylabel(r"dw-GNM Kirchhoff Fundamental $\nu_1$ (THz)", labelpad=8)
    ax1.set_xlim(0.01, lim_max)
    ax1.set_ylim(0.01, lim_max)
    ax1.legend(loc="upper left", fontsize=9.8, facecolor="#FFFFFF", edgecolor="#000000", framealpha=0.96)

    # =========================================================================
    # Panel B: Power-Law Scaling Curve \nu_1 vs \langle z \rangle / (R_g / R_0)
    # =========================================================================
    ax2 = axes[0, 1]
    ax2.grid(True)
    ax2.scatter(z_rg[is_enz], nu_calc[is_enz], s=26, alpha=0.68, color="#008080", label="Catalytic Enzymes (N = 922)")
    ax2.scatter(z_rg[~is_enz], nu_calc[~is_enz], s=38, alpha=0.90, color="#D62828", marker="^", label="Structural / Fibrillar (N = 118)")
    x_curve = np.linspace(float(np.min(z_rg)) * 0.95, float(np.max(z_rg)) * 1.03, 300)
    y_full = 0.42338 * (x_curve ** 1.1194)
    y_low = 0.4115 * (x_curve ** 1.1611)
    y_high = 0.4348 * (x_curve ** 1.0736)
    y_enz = 0.42951 * (x_curve ** 1.1570)
    y_str = 0.39597 * (x_curve ** 0.9037)
    ax2.plot(x_curve, y_full, color="#000000", lw=2.5, label=r"Global Fit: $\nu_1 = 0.423 (\langle z \rangle / (R_g/R_0))^{1.12}$ ($R^2 = 0.805$)")
    ax2.plot(x_curve, y_enz, color="#005F5F", lw=2.0, ls=":", label=r"Enzyme Regime: $b = 1.16$ ($R^2 = 0.799$)")
    ax2.plot(x_curve, y_str, color="#9B111E", lw=2.0, ls="-.", label=r"Structural Regime: $b = 0.90$ ($R^2 = 0.870$)")
    ax2.fill_between(x_curve, np.minimum(y_low, y_high), np.maximum(y_low, y_high), color="#64748B", alpha=0.22, label="95% Bootstrap CI Band")

    eq_text = (
        "Governing Scaling Law (N = 1,040):\n"
        r"$\nu_1\ (\mathrm{THz}) = 0.423 \left(\frac{\langle z \rangle}{R_g / R_0}\right)^{1.12}$" + "\n"
        r"$R^2 = 0.805 \quad|\quad r = 0.895 \quad|\quad p = 0.0002$"
    )
    ax2.text(
        0.965, 0.14, eq_text,
        transform=ax2.transAxes,
        ha="right", va="bottom",
        fontsize=10.6,
        fontweight="bold",
        color="#000000",
        bbox=dict(boxstyle="round,pad=0.55", facecolor="#F8FAFC", edgecolor="#000000", lw=1.4, alpha=0.98),
    )

    ax2.set_title("B. Power-Law Scaling & Functional Elasticity Regimes", pad=12)
    ax2.set_xlabel(r"Dimensionless 3D Contact Density $\langle z \rangle / (R_g / R_0)$ ($R_0 = 1.0\ \mathrm{\AA}$)", labelpad=7)
    ax2.set_ylabel(r"dw-GNM Kirchhoff Fundamental $\nu_1$ (THz)", labelpad=8)
    ax2.legend(loc="upper left", fontsize=9.5, facecolor="#FFFFFF", edgecolor="#000000", framealpha=0.96)

    # =========================================================================
    # Panel C: Variance Explained (R^2 %) Comparison
    # =========================================================================
    ax3 = axes[1, 0]
    ax3.grid(True, axis="y")
    categories = [
        "Octave Ctrl\n(mod-2)\n(N = 1,040)",
        "1D Length\n" + r"$1/\sqrt{N_{C\alpha}}$" + "\n(N = 115)",
        "1D Length\n" + r"$1/\sqrt{N_{C\alpha}}$" + "\n(N = 1,040)",
        "3D Radius\n" + r"$1/R_g$" + "\n(N = 1,040)",
        "3D Topology\n" + r"$\langle z \rangle / R_g$" + "\n(N = 1,040)",
        "3D Topology\n" + r"$\langle z \rangle / R_g$" + "\n(N = 115)",
    ]
    r2_values = [0.06, 3.2, 16.0, 66.3, 80.5, 81.2]
    r2_labels = ["0.06%", "3.2%", "16.0%", "66.3%", "80.5%", "81.2%"]
    bar_colors = ["#94A3B8", "#64748B", "#D97706", "#0284C7", "#008080", "#0F766E"]
    bars = ax3.bar(categories, r2_values, color=bar_colors, edgecolor="#000000", linewidth=1.0, width=0.52)
    for b_rect, val, lbl in zip(bars, r2_values, r2_labels):
        ax3.text(
            b_rect.get_x() + b_rect.get_width() / 2.0,
            val + 1.8,
            lbl,
            ha="center", va="bottom", fontsize=11.0, fontweight="bold", color="#000000"
        )
    ax3.set_ylim(0, 95)
    ax3.set_title(r"C. Variance Explained ($R^2$): 3D Topology vs. 1D Length", pad=12)
    ax3.set_ylabel(r"Variance in dw-GNM $\nu_1$ Explained ($R^2$ %)", labelpad=8)
    ax3.tick_params(axis="x", labelsize=10.2, pad=5)

    # =========================================================================
    # Panel D: Out-of-Sample Generalization & Kingdom Transfer
    # =========================================================================
    ax4 = axes[1, 1]
    ax4.grid(True, axis="x")
    protocols = [
        "Kingdom Transfer:\nMammal -> Non-Mammal (N = 911)",
        "Catalytic Enzymes\nSub-Cohort (N = 922)",
        "10-Fold Cross-Validation\nPooled Out-of-Sample (N = 1,040)",
        "Full Expanded Cohort\nPower-Law Fit (N = 1,040)",
        "80/20 Train/Test Split\nHeld-Out Test Set (N = 208)",
        "Primary Cohort\n10-Fold CV (N = 115)",
        "Kingdom Transfer:\nNon-Mammal -> Mammal (N = 129)",
        "Structural / Fibrillar\nSub-Cohort (N = 118)",
    ]
    oos_r2 = [78.4, 79.9, 80.4, 80.5, 80.5, 81.2, 86.1, 87.0]
    p_colors = ["#0284C7", "#008080", "#0D9488", "#0F766E", "#15803D", "#008080", "#1D4ED8", "#D62828"]
    y_pos = np.arange(len(protocols))
    hbars = ax4.barh(y_pos, oos_r2, color=p_colors, edgecolor="#000000", linewidth=0.9, height=0.62)
    ax4.set_yticks(y_pos)
    ax4.set_yticklabels(protocols, fontsize=10.2, fontweight="bold", color="#000000")
    for hb, val in zip(hbars, oos_r2):
        ax4.text(val + 1.0, hb.get_y() + hb.get_height() / 2.0, f"{val:.1f}%", va="center", fontsize=10.8, fontweight="bold", color="#000000")
    ax4.set_xlim(0, 99)
    ax4.set_title(r"D. Out-of-Sample Generalization ($R^2_{\mathrm{OOS}}$ %)", pad=12)
    ax4.set_xlabel(r"Out-of-Sample / Regime Variance Explained ($R^2_{\mathrm{OOS}}$ %)", labelpad=7)

    plt.tight_layout(rect=[0.012, 0.012, 0.988, 0.925], h_pad=2.4, w_pad=2.4)

    for target_dir in [LEGACY_DIR, ACTIVE_DIR]:
        if target_dir.exists():
            png_path = target_dir / "fig1_n1040_headline_scaling_law.png"
            pdf_path = target_dir / "fig1_n1040_headline_scaling_law.pdf"
            tif_path = target_dir / "fig1_n1040_headline_scaling_law.tif"
            plt.savefig(png_path, dpi=600, facecolor="#FFFFFF")
            plt.savefig(pdf_path, format="pdf", facecolor="#FFFFFF")
            plt.savefig(tif_path, dpi=600, facecolor="#FFFFFF")
            print(f"Saved White-Background 600-DPI PNG, Vector PDF, and 600-DPI TIFF to {target_dir}")

    plt.close(fig)


if __name__ == "__main__":
    main()
