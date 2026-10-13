#!/usr/bin/env python3
"""
07_figures.py - figures from results/. Nothing is hardcoded: every annotated statistic
is read from the JSON/CSV written upstream, which the previous version did not do for
the Trpc5 p value.

Fig 1  population definition and effector balance (snRNA)
Fig 2  Trpc5 depth control and the discreteness test   [new]
Fig 3  translatome, IP over Sup                         [new enrichment panel]
Fig 4  cross-population placement (Rupp)
"""
import os, json, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
import config as C
import vmhlib as V

plt.rcParams.update({"figure.dpi": 300, "savefig.dpi": 300, "font.size": 9,
                     "axes.titlesize": 10, "axes.titleweight": "bold",
                     "axes.spines.top": False, "axes.spines.right": False,
                     "font.family": "DejaVu Sans"})
R = C.RESULTS


def j(f): return json.load(open(os.path.join(R, f)))
def t(f):
    with open(os.path.join(R, f)) as fh:
        return list(csv.DictReader(fh))
def it(g): return f"$\\it{{{g}}}$"
def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(C.FIGURES, f"{name}.{ext}"), bbox_inches="tight")


def fig1():
    det, cs = j("detection.json"), t("cluster_scores.csv")
    pc = j("percell_stats.json")
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))

    g = det["naive_gate"]
    n_lepr = det["n_lepr_labelled"]
    genes = ["Nr5a1", "Lepr", "Slc17a6"]
    x = np.arange(len(genes)); w = .38
    v1 = [det["vmh_detect"][k] for k in genes]; v2 = [det["lepr_detect"][k] for k in genes]
    ax[0].bar(x - w/2, v1, w, label="all VMH neurons", color=C.COLORS["NEUTRAL"])
    ax[0].bar(x + w/2, v2, w, label="LepR clusters", color=C.COLORS["LEPR"])
    ax[0].set_xticks(x); ax[0].set_xticklabels([it(k) for k in genes])
    ax[0].set_ylabel("% nuclei with >=1 UMI"); ax[0].set_ylim(0, 100)
    ax[0].set_title("A  Dropout in the defining transcripts"); ax[0].legend(loc="upper left")

    ax[1].bar(["cluster-defined", "Nr5a1+ and Lepr+"], [n_lepr, g["lepr_clusters"]],
              color=[C.COLORS["LEPR"], "#999999"], edgecolor="white")
    for i, v in enumerate([n_lepr, g["lepr_clusters"]]):
        ax[1].text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=8)
    ax[1].set_ylabel("LepR nuclei retained")
    ax[1].set_title(f"B  A per-cell gate keeps {g['lepr_pct']}%")

    rr = sorted(cs, key=lambda r: float(r["HminusD"]))
    y = np.arange(len(rr))
    cols = [C.COLORS["LEPR"] if r["is_lepr"] == "True" else C.COLORS["OTHER"] for r in rr]
    ax[2].barh(y, [float(r["HminusD"]) for r in rr], color=cols, edgecolor="white", lw=.3)
    ax[2].set_yticks(y); ax[2].set_yticklabels([r["cluster"] for r in rr], fontsize=5.5)
    ax[2].axvline(0, color="k", lw=.6)
    ax[2].set_xlabel("hyperpolarizing - depolarizing (log-CP10k)")
    k = pc["HminusD"]
    ax[2].set_title(f"C  Effector axis (per-cell eps2 = {k['eps2_mean']}, "
                    f"{k['eps2_min']} to {k['eps2_max']})")
    ax[2].legend(handles=[Patch(color=C.COLORS["LEPR"], label="LepR"),
                          Patch(color=C.COLORS["OTHER"], label="other VMH")],
                 loc="lower right", fontsize=7)
    fig.tight_layout(); save(fig, "Fig1_population_effector"); return fig


def fig2():
    d, disc = j("trpc5_depth.json"), j("discreteness.json")
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))

    raw = d["downsampling"]["pct_detected_raw"]
    ds = d["downsampling"]["pct_detected_downsampled"]
    ks = sorted(raw, key=lambda k: -raw[k]); x = np.arange(len(ks)); w = .38
    ax[0].bar(x - w/2, [raw[k] for k in ks], w, label="observed", color=C.COLORS["DEPOL"])
    ax[0].bar(x + w/2, [ds[k] for k in ks], w, label=f"downsampled to "
              f"{d['downsampling']['target_umi']} UMI", color="#B0B0B0")
    ax[0].set_xticks(x); ax[0].set_xticklabels(ks, rotation=90, fontsize=6)
    ax[0].set_ylabel(f"% nuclei {it('Trpc5')}+")
    lg = d["logistic_cluster_given_depth"]
    ax[0].set_title(f"A  {it('Trpc5')} with and without depth control\n"
                    f"LR(cluster|depth) = {lg['LR']}, p = {lg['p']:.1e}", fontsize=8.5)
    ax[0].legend(fontsize=6.5)

    null = [float(r["lr"]) for r in t("trpc5_depth_null.csv")]
    ax[1].hist(null, bins=30, color="#B0B0B0", edgecolor="white")
    ax[1].axvline(lg["LR"], color=C.COLORS["DEPOL"], lw=1.5)
    ax[1].text(lg["LR"], ax[1].get_ylim()[1] * .9, f" {it('Trpc5')}", color=C.COLORS["DEPOL"],
               fontsize=7, va="top")
    ax[1].set_xlabel("LR statistic for cluster, given depth")
    ax[1].set_ylabel("depth-matched control genes")
    ax[1].set_title(f"B  Empirical null\n{it('Trpc5')} exceeds "
                    f"{d['empirical_null']['trpc5_percentile']:.0f}% of controls", fontsize=8.5)

    lev = disc["power_simulation"]["levels"]
    ax[2].plot([l["effect_sd_units"] for l in lev],
               [l["pct_sims_bic_picks_k_ge_2"] for l in lev], "o-", color=C.COLORS["HYPER"])
    ax[2].axhline(80, color="0.6", ls="--", lw=.7)
    ax[2].set_xlabel("simulated split (SD units)")
    ax[2].set_ylabel("% simulations where BIC selects k >= 2")
    ax[2].set_ylim(0, 105)
    ax[2].set_title(f"C  Power to see a discrete split\nobserved best k = "
                    f"{disc['observed']['best_k_by_bic']}", fontsize=8.5)
    fig.tight_layout(); save(fig, "Fig2_depth_discreteness"); return fig


def fig3():
    cpm, enr = t("trap_cpm.csv"), t("trap_enrichment.csv")
    keep = [r["gene"] for r in cpm if r["gene"] in C.KEY_CHANNELS]
    cd = {r["gene"]: r for r in cpm}; ed = {r["gene"]: r for r in enr}
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
    x = np.arange(len(keep)); w = .38
    ax[0].bar(x - w/2, [float(cd[g]["SF1_Bead_mean"]) for g in keep], w,
              color=C.COLORS["TRAP"], label="SF1 Bead (IP)")
    ax[0].bar(x + w/2, [float(cd[g]["SF1_Sup_mean"]) for g in keep], w,
              color="#B8D8C4", edgecolor=C.COLORS["TRAP"], label="SF1 Sup")
    ax[0].set_xticks(x); ax[0].set_xticklabels([it(g) for g in keep], rotation=40, ha="right",
                                               fontsize=7)
    ax[0].set_ylabel("CPM"); ax[0].legend(fontsize=7)
    ax[0].set_title("A  TRAP fractions, SF1 (n = 3 per fraction)")

    vals = [float(ed[g]["SF1_log2_IP_over_Sup"]) for g in keep]
    cols = [C.COLORS["HYPER"] if g in C.PANELS["HYPER"] else C.COLORS["DEPOL"] for g in keep]
    ax[1].barh(np.arange(len(keep)), vals, color=cols, edgecolor="white")
    ax[1].axvline(0, color="k", lw=.6)
    ax[1].set_yticks(np.arange(len(keep))); ax[1].set_yticklabels([it(g) for g in keep],
                                                                 fontsize=7)
    ax[1].set_xlabel("log2 IP / Sup")
    ax[1].set_title("B  Enrichment, not bare IP abundance")
    ax[1].legend(handles=[Patch(color=C.COLORS["HYPER"], label="hyperpolarizing"),
                          Patch(color=C.COLORS["DEPOL"], label="depolarizing")], fontsize=7)
    fig.tight_layout(); save(fig, "Fig3_TRAP_enrichment"); return fig


def fig4():
    sc, xp = t("crosspop_scores.csv"), j("crosspop.json")
    pops = [r["population"] for r in sc]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.0))
    for i, (lab, key) in enumerate((("H - D axis", "HminusD"), (f"{it('Trpc5')}", "Trpc5"))):
        m = [float(r[f"{key}_mean"]) for r in sc]
        lo = [float(r[f"{key}_mean"]) - float(r[f"{key}_ci_lo"]) for r in sc]
        hi = [float(r[f"{key}_ci_hi"]) - float(r[f"{key}_mean"]) for r in sc]
        ax[i].bar(pops, m, yerr=[lo, hi], capsize=4,
                  color=[C.COLORS["LEPR"] if p == "VMH_SF1" else "#7F7F7F" for p in pops])
        ax[i].axhline(0, color="k", lw=.6)
        ax[i].set_ylabel(lab)
        ax[i].set_title(("A  " if i == 0 else "B  ") + lab +
                        "  (bootstrap 95% CI)", fontsize=9)
    p = xp["contrasts"]["POMC_vs_AgRP_Trpc5"]["perm_p_seed1"]
    ax[1].text(.02, .95, f"POMC > AgRP, permutation p = {p:.4f}", transform=ax[1].transAxes,
               fontsize=7, va="top")
    fig.tight_layout(); save(fig, "Fig4_crosspopulation"); return fig


def main():
    figs = [fig1(), fig2(), fig3(), fig4()]
    with PdfPages(os.path.join(C.FIGURES, "VMH_leptin_figures.pdf")) as pdf:
        for f in figs:
            pdf.savefig(f, bbox_inches="tight")
    print("Figures in", C.FIGURES)
    for f in sorted(os.listdir(C.FIGURES)):
        print("  ", f)


if __name__ == "__main__":
    main()
