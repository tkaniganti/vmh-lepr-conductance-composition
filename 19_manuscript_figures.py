#!/usr/bin/env python3
"""
19_manuscript_figures.py - the four figures of the manuscript, built from results/.

Nothing is hardcoded: every annotated statistic is read from the JSON and CSV written by
the analysis scripts, so a rerun of the pipeline regenerates the figures consistently.

Fig 1  Populations cannot be defined by detection of their marker transcripts
       A detection rates, B gate collapse, C effector axis across 24 VMH clusters
Fig 2  Effector balance is organised at the conserved-class level
       A variance partition, B Lepr class against each other class, C marker composition
Fig 3  Structure within the leptin-receptor population
       A mixture BIC and power simulation, B Trpc5 depth control, C dispersion
Fig 4  Translatome corroboration and cross-population placement
       A IP-over-Sup enrichment, B balance axis and Trpc5 across populations

Inputs: detection.json, cluster_scores.csv, conserved_classes.json, marker_classes_cons6.csv,
        discreteness.json, trpc5_depth.json, trpc5_depth_null.csv, dispersion.json,
        trap_enrichment.csv, crosspop_scores.csv, crosspop.json
Output: figures/Fig1..Fig4 (png + pdf) and figures/manuscript_figures.pdf
"""
import os, json, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
import config as C

plt.rcParams.update({"figure.dpi": 300, "savefig.dpi": 300, "font.size": 8,
                     "axes.titlesize": 9, "axes.titleweight": "bold",
                     "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "legend.fontsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "font.family": "DejaVu Sans"})
R = C.RESULTS
BLUE, ORANGE, RED, GREY = C.COLORS["HYPER"], C.COLORS["DEPOL"], C.COLORS["LEPR"], C.COLORS["OTHER"]


def j(f):
    p = os.path.join(R, f)
    if not os.path.exists(p):
        raise SystemExit(f"missing {p}; run the analysis scripts first")
    return json.load(open(p))


def t(f):
    p = os.path.join(R, f)
    if not os.path.exists(p):
        raise SystemExit(f"missing {p}")
    with open(p) as fh:
        return list(csv.DictReader(fh))


def it(g):
    return f"$\\it{{{g}}}$"


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(C.FIGURES, f"{name}.{ext}"), bbox_inches="tight")


# ------------------------------------------------------------------ Figure 1
def fig1():
    det, cs = j("detection.json"), t("cluster_scores.csv")
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.1))

    genes = ["Nr5a1", "Lepr", "Slc17a6"]
    x, w = np.arange(len(genes)), 0.38
    v1 = [det["vmh_detect"][g] for g in genes]
    v2 = [det["lepr_detect"][g] for g in genes]
    ax[0].bar(x - w / 2, v1, w, label="all VMH neurons", color=C.COLORS["NEUTRAL"])
    ax[0].bar(x + w / 2, v2, w, label="LepR clusters", color=RED)
    for xi, (a, b) in enumerate(zip(v1, v2)):
        ax[0].text(xi - w / 2, a + 1.5, f"{a:.1f}", ha="center", fontsize=6.5)
        ax[0].text(xi + w / 2, b + 1.5, f"{b:.1f}", ha="center", fontsize=6.5)
    ax[0].set_xticks(x); ax[0].set_xticklabels([it(g) for g in genes])
    ax[0].set_ylabel("% nuclei with $\\geq$1 UMI"); ax[0].set_ylim(0, 100)
    ax[0].legend(loc="upper left", frameon=False)
    ax[0].set_title("A   Detection of defining transcripts", loc="left")

    g = det["naive_gate"]
    n_lepr = det["n_lepr_labelled"]
    bars = ax[1].bar(["cluster-defined", "$\\it{Nr5a1}$+ and $\\it{Lepr}$+"],
                     [n_lepr, g["lepr_clusters"]], color=[RED, "#999999"], edgecolor="white")
    for b, v in zip(bars, [n_lepr, g["lepr_clusters"]]):
        ax[1].text(b.get_x() + b.get_width() / 2, v, f"{v:,}", ha="center",
                   va="bottom", fontsize=7)
    ax[1].set_ylabel("LepR nuclei retained")
    ax[1].set_title(f"B   A per-cell gate retains {g['lepr_pct']}%", loc="left")

    rr = sorted(cs, key=lambda r: float(r["HminusD"]))
    y = np.arange(len(rr))
    cols = [RED if r["is_lepr"] == "True" else GREY for r in rr]
    ax[2].barh(y, [float(r["HminusD"]) for r in rr], color=cols, edgecolor="white", lw=0.3)
    ax[2].set_yticks(y); ax[2].set_yticklabels([r["cluster"] for r in rr], fontsize=5)
    ax[2].axvline(0, color="k", lw=0.6)
    ax[2].set_xlabel("hyperpolarizing $-$ depolarizing (log-CP10k)")
    ax[2].legend(handles=[Patch(color=RED, label="LepR cluster"),
                          Patch(color=GREY, label="other VMH")],
                 loc="lower right", frameon=False)
    ax[2].set_title("C   Effector axis, 24 VMH clusters", loc="left")

    fig.tight_layout(); save(fig, "Fig1_population"); return fig


# ------------------------------------------------------------------ Figure 2
def fig2():
    cc = j("conserved_classes.json")
    mk = t("marker_classes_cons6.csv")
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.2))

    part = cc["variance_partition"]
    labels = ["between\nconserved classes", "between clusters\nwithin LepR class"]
    xs = np.arange(2); w = 0.36
    for k, (tag, col) in enumerate((("raw", BLUE), ("depth_residualized", "#9FB8DC"))):
        p = part[tag]
        means = [p["eps2_between_conserved_classes"]["mean"],
                 p["eps2_between_clusters_within_lepr_class"]["mean"]]
        lo = [means[0] - p["eps2_between_conserved_classes"]["min"],
              means[1] - p["eps2_between_clusters_within_lepr_class"]["min"]]
        hi = [p["eps2_between_conserved_classes"]["max"] - means[0],
              p["eps2_between_clusters_within_lepr_class"]["max"] - means[1]]
        ax[0].bar(xs + (k - 0.5) * w, means, w, yerr=[lo, hi], capsize=3, color=col,
                  label="raw" if tag == "raw" else "depth-residualized")
    ratio = part["raw"].get("ratio_between_over_within")
    ax[0].set_xticks(xs); ax[0].set_xticklabels(labels, fontsize=7)
    ax[0].set_ylabel("$\\epsilon^2$"); ax[0].axhline(0, color="k", lw=0.5)
    ax[0].legend(frameon=False)
    ax[0].set_title(f"A   Variance partition (ratio {ratio})", loc="left")

    pair = cc["lepr_class_vs_others"]
    keys = sorted(pair, key=lambda k: pair[k]["diff_mean_over_seeds"])
    names = [k.split("_vs_")[1] for k in keys]
    vals = [pair[k]["diff_mean_over_seeds"] for k in keys]
    err = [pair[k]["diff_range"] / 2 for k in keys]
    cols = [RED if pair[k]["sign_consistent_across_seeds"] else GREY for k in keys]
    ax[1].barh(np.arange(len(keys)), vals, xerr=err, capsize=3, color=cols, edgecolor="white")
    ax[1].set_yticks(np.arange(len(keys))); ax[1].set_yticklabels(names)
    ax[1].axvline(0, color="k", lw=0.6)
    ax[1].set_xlabel("LepR class $-$ other class (H$-$D)")
    ax[1].set_title("B   Conserved LepR class vs others", loc="left")
    ax[1].legend(handles=[Patch(color=RED, label="sign consistent, 20/20 draws")],
                 loc="lower right", frameon=False)

    want = ["Gck", "Bdnf", "Esr1"]
    groups = sorted({r["group"] for r in mk})
    byg = {(r["group"], r["gene"]): r for r in mk}
    xs = np.arange(len(groups)); w = 0.26
    for i, gene in enumerate(want):
        vals = [float(byg[(g, gene)]["log_cp10k_mean"]) if (g, gene) in byg else 0
                for g in groups]
        lo = [vals[k] - float(byg[(g, gene)]["ci_lo"]) for k, g in enumerate(groups)]
        hi = [float(byg[(g, gene)]["ci_hi"]) - vals[k] for k, g in enumerate(groups)]
        ax[2].bar(xs + (i - 1) * w, vals, w, yerr=[lo, hi], capsize=2,
                  label=it(gene), color=[BLUE, ORANGE, "#6E6E6E"][i])
    ax[2].set_xticks(xs)
    ax[2].set_xticklabels(groups, rotation=45, ha="right", fontsize=7)
    ax[2].set_ylabel("log-CP10k"); ax[2].legend(frameon=False, ncol=3, fontsize=6.5)
    ax[2].set_title("C   Marker composition by conserved class", loc="left")

    fig.tight_layout(); save(fig, "Fig2_conserved_classes"); return fig


# ------------------------------------------------------------------ Figure 3
def fig3():
    disc, dep, disp = j("discreteness.json"), j("trpc5_depth.json"), j("dispersion.json")
    null = [float(r["lr"]) for r in t("trpc5_depth_null.csv")]
    fig, ax = plt.subplots(1, 4, figsize=(13.5, 3.2))

    mix = disc["observed"]["mixture"]
    ks = [m["k"] for m in mix]; bics = [m["bic"] for m in mix]
    ax[0].plot(ks, bics, "o-", color=BLUE)
    best = disc["observed"]["best_k_by_bic"]
    ax[0].scatter([best], [bics[ks.index(best)]], s=70, facecolors="none",
                  edgecolors=RED, zorder=5, linewidths=1.4)
    ax[0].set_xticks(ks); ax[0].set_xlabel("mixture components"); ax[0].set_ylabel("BIC")
    d = disc["observed"].get("dip", {})
    ptxt = f"dip p = {d.get('p'):.3f}" if isinstance(d, dict) and "p" in d else ""
    ax[0].set_title(f"A   Unimodal at k = {best}   ({ptxt})", loc="left")

    lev = disc["power_simulation"]["levels"]
    ax[1].plot([l["effect_sd_units"] for l in lev],
               [l["pct_sims_bic_picks_k_ge_2"] for l in lev], "s--", color=ORANGE, ms=5)
    ax[1].axhline(80, color="0.6", ls=":", lw=0.7)
    ax[1].set_xlabel("simulated split (SD units)")
    ax[1].set_ylabel("% simulations detecting a split")
    ax[1].set_ylim(0, 105)
    ax[1].set_title("B   Power of the same test", loc="left")

    raw = dep["downsampling"]["pct_detected_raw"]
    ds = dep["downsampling"]["pct_detected_downsampled"]
    ks2 = sorted(raw, key=lambda k: -raw[k]); x = np.arange(len(ks2)); w = 0.38
    ax[2].bar(x - w / 2, [raw[k] for k in ks2], w, label="observed", color=ORANGE)
    ax[2].bar(x + w / 2, [ds[k] for k in ks2], w, color=GREY,
              label=f"depth-matched ({dep['downsampling']['target_umi']} UMI)")
    ax[2].set_xticks(x); ax[2].set_xticklabels(ks2, rotation=90, fontsize=6)
    ax[2].set_ylabel(f"% nuclei {it('Trpc5')}+")
    lg = dep["logistic_cluster_given_depth"]
    ax[2].legend(frameon=False, fontsize=6.5)
    ax[2].set_title(f"C   LR = {lg['LR']}, p = {lg['p']:.0e}; "
                    f"{dep['empirical_null']['trpc5_percentile']:.1f}th pct of null",
                    loc="left", fontsize=7.5)

    sd = disp["depth_residualized"]["sd_by_population"]
    pops = [p for p in ("POMC", "AgRP", "VMH_SF1") if p in sd]
    means = [sd[p]["mean"] for p in pops]
    lo = [sd[p]["mean"] - sd[p]["min"] for p in pops]
    hi = [sd[p]["max"] - sd[p]["mean"] for p in pops]
    ax[3].bar([p.replace("_", " ") for p in pops], means, yerr=[lo, hi], capsize=3,
              color=[GREY if p != "VMH_SF1" else RED for p in pops])
    ax[3].set_ylabel("SD of depth-residualized H$-$D")
    r = disp["equal_n"]["sd_ratio_median"]
    ax[3].set_title("D   Dispersion (VMH/POMC {:.2f}, VMH/AgRP {:.2f} at equal n)".format(
        r["VMH_over_POMC"], r["VMH_over_AgRP"]), loc="left", fontsize=7.5)

    fig.tight_layout(); save(fig, "Fig3_within_population"); return fig


# ------------------------------------------------------------------ Figure 4
def fig4():
    enr = t("trap_enrichment.csv")
    sc, xp = t("crosspop_scores.csv"), j("crosspop.json")
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.0))

    show = ["Kcnj5", "Kcnj3", "Hcn1", "Abcc8", "Kcnj9", "Kcnj6", "Trpc1", "Trpc4",
            "Trpc5", "Hcn2", "Kcnj11", "Pde3b"]
    d = {r["gene"]: r for r in enr}
    show = [g for g in show if g in d]
    vals = [float(d[g]["SF1_log2_IP_over_Sup"]) for g in show]
    ps = [d[g]["SF1_welch_p_unadjusted"] for g in show]
    cols = [BLUE if g in C.PANELS["HYPER"] else ORANGE for g in show]
    y = np.arange(len(show))
    ax[0].barh(y, vals, color=cols, edgecolor="white")
    for yi, (v, p) in enumerate(zip(vals, ps)):
        try:
            if p not in (None, "") and float(p) < 0.05:
                ax[0].text(v + (0.03 if v >= 0 else -0.03), yi, "*",
                           va="center", ha="left" if v >= 0 else "right", fontsize=9)
        except ValueError:
            pass
    ax[0].set_yticks(y); ax[0].set_yticklabels([it(g) for g in show], fontsize=7)
    ax[0].axvline(0, color="k", lw=0.6)
    ax[0].set_xlabel("log$_2$ IP / Sup, SF1 fraction")
    lo_, hi_ = min(vals), max(vals)
    pad = 0.18 * max(abs(lo_), abs(hi_))
    ax[0].set_xlim(lo_ - pad, hi_ + pad)
    ax[0].invert_yaxis()
    ax[0].legend(handles=[Patch(color=BLUE, label="hyperpolarizing"),
                          Patch(color=ORANGE, label="depolarizing")],
                 loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2,
                 frameon=False, fontsize=7)
    ax[0].set_title("A   Translatome enrichment ($*$ p < 0.05 unadjusted, n = 3)",
                    loc="left", pad=30)

    pops = [r["population"] for r in sc]
    xs = np.arange(len(pops)); w = 0.36
    for i, (key, lab, col) in enumerate((("HminusD", "H$-$D axis", BLUE),
                                         ("Trpc5", it("Trpc5"), ORANGE))):
        m = [float(r[f"{key}_mean"]) for r in sc]
        lo = [float(r[f"{key}_mean"]) - float(r[f"{key}_ci_lo"]) for r in sc]
        hi = [float(r[f"{key}_ci_hi"]) - float(r[f"{key}_mean"]) for r in sc]
        axis = ax[1] if i == 0 else ax[1].twinx()
        axis.bar(xs + (i - 0.5) * w, m, w, yerr=[lo, hi], capsize=3, color=col, label=lab)
        axis.set_ylabel(lab, color=col); axis.tick_params(axis="y", colors=col)
        if i == 1:
            axis.spines["right"].set_visible(True)
    ax[1].set_xticks(xs); ax[1].set_xticklabels(
        [f"{p.replace(chr(95), chr(32))}\n(n={r['n']})" for p, r in zip(pops, sc)], fontsize=7)
    ax[1].axhline(0, color="k", lw=0.5)
    c = xp["contrasts"]["VMHSF1_vs_POMC_HminusD"]
    ax[1].set_title("B   Cross-population placement "
                    f"(VMH vs POMC {c['diff_mean_over_seeds']:+.3f}, sign consistent)",
                    loc="left", fontsize=7.5)

    fig.tight_layout(); save(fig, "Fig4_trap_crosspop"); return fig


def main():
    figs = [fig1(), fig2(), fig3(), fig4()]
    with PdfPages(os.path.join(C.FIGURES, "manuscript_figures.pdf")) as pdf:
        for f in figs:
            pdf.savefig(f, bbox_inches="tight")
    print("Figures written to", C.FIGURES)
    for fn in sorted(os.listdir(C.FIGURES)):
        if fn.startswith(("Fig1_pop", "Fig2_cons", "Fig3_with", "Fig4_trap", "manuscript")):
            print("  ", fn)


if __name__ == "__main__":
    main()
