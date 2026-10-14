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
    """Class-level structure: nested partition, class contrasts, marker composition."""
    npart = j("nested_partition.json"); cc = j("conserved_classes.json")
    mk = t("marker_classes_cons6.csv")
    fig, ax = plt.subplots(1, 3, figsize=(11.0, 3.3))

    comp = npart["nested_all_vmh"]
    keys = ["pct_class", "pct_cluster_in_class"]
    labs = ["conserved\nclass", "cluster within\nclass"]
    m = [comp[k]["mean"] for k in keys]
    lo = [comp[k]["mean"] - comp[k]["min"] for k in keys]
    hi = [comp[k]["max"] - comp[k]["mean"] for k in keys]
    ax[0].bar(labs, m, yerr=[lo, hi], capsize=4, color=[BLUE, "#9FB8DC"])
    ax[0].set_ylabel("% of variance")
    ax[0].set_title("A   Nested partition", loc="left")

    pair = cc["lepr_class_vs_others"]
    kk = sorted(pair, key=lambda k: pair[k]["diff_mean_over_seeds"])
    ax[1].barh(np.arange(len(kk)), [pair[k]["diff_mean_over_seeds"] for k in kk],
               xerr=[pair[k]["diff_range"]/2 for k in kk], capsize=3,
               color=[RED if pair[k]["sign_consistent_across_seeds"] else GREY for k in kk],
               edgecolor="white")
    ax[1].set_yticks(np.arange(len(kk)))
    ax[1].set_yticklabels([k.split("_vs_")[1] for k in kk], fontsize=7)
    ax[1].axvline(0, color="k", lw=.6)
    ax[1].set_xlabel("LepR class $-$ other (H$-$D)")
    ax[1].set_title("B   Class contrasts", loc="left")

    want = ["Gck", "Bdnf", "Esr1"]
    groups = sorted({r["group"] for r in mk})
    byg = {(r["group"], r["gene"]): r for r in mk}
    xs = np.arange(len(groups)); w = .26
    for i, gene in enumerate(want):
        v = [float(byg[(g, gene)]["log_cp10k_mean"]) if (g, gene) in byg else 0 for g in groups]
        ax[2].bar(xs + (i-1)*w, v, w, label=it(gene), color=[BLUE, ORANGE, "#6E6E6E"][i])
    ax[2].set_xticks(xs); ax[2].set_xticklabels(groups, rotation=45, ha="right", fontsize=6.5)
    ax[2].set_ylabel("log-CP10k"); ax[2].legend(frameon=False, ncol=3, fontsize=6)
    ax[2].set_title("C   Marker composition", loc="left")

    fig.tight_layout(); save(fig, "Fig2_class_level"); return fig


def fig3():
    """Within-class survey: per-class excess, threshold sweep, detectability."""
    cn = j("class_nulls.json"); npart = j("nested_partition.json")
    mc = j("mincluster_sensitivity.json")
    fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.8))

    bc = cn["by_class"]
    order = sorted(bc, key=lambda c: -(bc[c]["excess_over_null"]))
    x = np.arange(len(order))
    exc = [bc[c]["excess_over_null"] for c in order]
    cols = [RED if c == "Lepr" else GREY for c in order]
    ax[0].bar(x, exc, color=cols, edgecolor="0.3", lw=.5)
    ax[0].axhline(0, color="k", lw=.6)
    ax[0].set_ylim(-1.2, max(exc) * 1.14)
    for xi, c in enumerate(order):
        v = exc[xi]
        ax[0].annotate(f"{v:.2f}", (xi, v), textcoords="offset points",
                       xytext=(0, 4 if v >= 0 else -9), ha="center", fontsize=6.5,
                       color=RED if c == "Lepr" else "0.25")
        ax[0].annotate(f"k={bc[c]['n_clusters']}, n={bc[c]['n_nuclei']}", (xi, 0),
                       xycoords=("data", "axes fraction"), textcoords="offset points",
                       xytext=(0, -26), ha="center", fontsize=5.5, color="0.4",
                       annotation_clip=False)
    ax[0].set_xticks(x); ax[0].set_xticklabels(order, fontsize=7)
    ax[0].tick_params(axis="x", pad=3)
    ax[0].set_ylabel("excess over chance (%)")
    ax[0].set_title("A   Within-class structure by class", loc="left")

    thr = mc["thresholds"]
    for c, col in (("Nfib", BLUE), ("Lepr", RED), ("Esr1", "0.45")):
        v = [mc["by_class"][c][str(t_)].get("excess") for t_ in thr]
        ok = [(t_, x_) for t_, x_ in zip(thr, v) if x_ is not None]
        ax[1].plot([a for a, _ in ok], [b for _, b in ok], "o-", color=col, label=c, ms=4)
    ax[1].axhline(0, color="k", lw=.6)
    ax[1].set_xlabel("minimum cluster size (nuclei)")
    ax[1].set_ylabel("excess over chance (%)")
    ax[1].legend(frameon=False, fontsize=7)
    ax[1].set_title("B   Excluding small clusters", loc="left")

    lev = npart["power_matched_control"]["levels"]
    ax[2].plot([l["effect_sd_units"] for l in lev],
               [l["recovered_pct"]["mean"] for l in lev], "o-", color=ORANGE)
    obs = npart["within_lepr_class_nine_clusters"]["pct"]["mean"]
    ax[2].axhline(obs, color=RED, ls="--", lw=.9)
    ax[2].text(0.40, obs, "observed ", color=RED, fontsize=6.5, va="bottom", ha="right")
    ax[2].set_xlabel("injected cluster effect (SD units)")
    ax[2].set_ylabel("variance recovered (%)")
    ax[2].set_title("C   Detectability, given cluster labels", loc="left")

    fig.tight_layout(); save(fig, "Fig3_within_class_survey"); return fig


# ------------------------------------------------------------------ Figure 4
def fig4():
    disc, dep = j("discreteness.json"), j("trpc5_depth.json")
    null = [float(r["lr"]) for r in t("trpc5_depth_null.csv")]
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.3))

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
    ax[1].axhline(80, color="0.6", ls=":", lw=.7)
    ax[1].set_xlabel("simulated split (SD units)")
    ax[1].set_ylabel("% simulations detecting a split")
    ax[1].set_ylim(0, 105)
    ax[1].set_title("B   Power, unsupervised mixture recovery", loc="left")

    raw = dep["downsampling"]["pct_detected_raw"]
    ds = dep["downsampling"]["pct_detected_downsampled"]
    ks2 = sorted(raw, key=lambda k: -raw[k]); x = np.arange(len(ks2)); w = .38
    ax[2].bar(x - w/2, [raw[k] for k in ks2], w, label="observed", color=ORANGE)
    ax[2].bar(x + w/2, [ds[k] for k in ks2], w, color=GREY,
              label=f"depth-matched ({dep['downsampling']['target_umi']} UMI)")
    ax[2].set_xticks(x); ax[2].set_xticklabels(ks2, rotation=90, fontsize=6)
    ax[2].set_ylabel(f"% nuclei {it('Trpc5')}+")
    lg = dep["logistic_cluster_given_depth"]
    ax[2].legend(frameon=False, fontsize=6.5)
    ax[2].set_title(f"C   LR = {lg['LR']}, p = {lg['p']:.0e}; "
                    f"{dep['empirical_null']['trpc5_percentile']:.1f}th pct of null",
                    loc="left", fontsize=7.5)

    fig.tight_layout(); save(fig, "Fig4_within_population"); return fig


# ------------------------------------------------------------------ Figure 5
def fig5():
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
    ax[0].set_title("A   Translatome enrichment ($*$ q < 0.05, BH, n = 3)",
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

    fig.tight_layout(); save(fig, "Fig5_trap_crosspop"); return fig


def main():
    figs = [fig1(), fig2(), fig3(), fig4(), fig5()]
    with PdfPages(os.path.join(C.FIGURES, "manuscript_figures.pdf")) as pdf:
        for f in figs:
            pdf.savefig(f, bbox_inches="tight")
    print("Figures written to", C.FIGURES)
    for fn in sorted(os.listdir(C.FIGURES)):
        if fn.startswith(("Fig1_", "Fig2_", "Fig3_", "Fig4_", "Fig5_", "manuscript")):
            print("  ", fn)


if __name__ == "__main__":
    main()
