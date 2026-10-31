#!/usr/bin/env python3
"""
19_manuscript_figures.py - the five figures of the manuscript, built from results/.

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

Inputs: detection.json, cluster_axis_corrected.json, power_grid.json, conserved_classes.json, marker_classes_cons6.csv,
        discreteness.json, trpc5_depth.json, trpc5_depth_null.csv, dispersion.json,
        trap_enrichment.csv, crosspop_scores.csv, crosspop.json
Output: figures/Fig1..Fig5 (png + pdf), figures/manuscript_figures.pdf and
        figures/checksums.txt, which records the md5 of every file written here.
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


WRITTEN = []


def save(fig, name):
    WRITTEN.append(name)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(C.FIGURES, f"{name}.{ext}"), bbox_inches="tight")


# ------------------------------------------------------------------ Figure 1
def fig1():
    det = j("detection.json")
    mk = t("marker_classes_cons6.csv")
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
    ax[0].set_ylabel("Nuclei with $\\geq$1 UMI (%)"); ax[0].set_ylim(0, 100)
    ax[0].legend(loc="upper left", frameon=False)
    ax[0].set_title("A   Detection of defining transcripts", loc="left")

    g = det["naive_gate"]
    # n_lepr_labelled is 2,914, the metadata count BEFORE the barcode filter; the
    # manuscript and the gate denominator both use 2,913 (n_lepr_with_counts). The
    # panel was printing 2,914 beside a gated bar of 316 of 2,913, so the two bars had
    # different denominators. Found in the rendered PDF, not by any check: a number
    # drawn inside a panel is compared to nothing. The manifest checks declarations
    # against results/, 44 checks declarations against the text, 39 checks text against
    # declarations. Panel annotations sit outside all three.
    n_lepr = det["n_lepr_with_counts"]
    bars = ax[1].bar(["cluster-defined", "$\\it{Nr5a1}$+ and $\\it{Lepr}$+"],
                     [n_lepr, g["lepr_clusters"]], color=[RED, "#999999"], edgecolor="white")
    for b, v in zip(bars, [n_lepr, g["lepr_clusters"]]):
        ax[1].text(b.get_x() + b.get_width() / 2, v, f"{v:,}", ha="center",
                   va="bottom", fontsize=7)
    ax[1].set_ylabel("LepR nuclei retained (n)")
    ax[1].set_title("B   Nuclei retained by a per-cell gate", loc="left")

    want = ["Gck", "Bdnf", "Esr1"]
    groups = sorted({r["group"] for r in mk})
    byg = {(r["group"], r["gene"]): r for r in mk}
    xs = np.arange(len(groups)); w = .26
    for i, gene in enumerate(want):
        v = [float(byg[(g, gene)]["log_cp10k_mean"]) if (g, gene) in byg else 0
             for g in groups]
        ax[2].bar(xs + (i - 1) * w, v, w, label=it(gene),
                  color=[BLUE, ORANGE, "#6E6E6E"][i])
    ax[2].set_xticks(xs); ax[2].set_xticklabels(groups, rotation=45, ha="right",
                                                fontsize=6.5)
    ax[2].set_ylabel("Expression (log-CP10k)")
    ax[2].legend(frameon=False, ncol=3, fontsize=6)
    ax[2].set_title("C   Marker expression by class", loc="left")

    fig.tight_layout(); save(fig, "Fig1_population"); return fig


# ------------------------------------------------------------------ Figure 2
def fig2():
    """Class-level structure: cluster ordering on the corrected axis, class contrasts."""
    cax = j("cluster_axis_corrected.json"); cc = j("conserved_classes.json")
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 4.2))

    bc = cax["by_cluster"]
    rr = sorted(bc, key=lambda c: bc[c]["HminusD"]["mean"])
    yv = np.arange(len(rr))
    mean = [bc[c]["HminusD"]["mean"] for c in rr]
    lo = [bc[c]["HminusD"]["mean"] - bc[c]["HminusD"]["min"] for c in rr]
    hi = [bc[c]["HminusD"]["max"] - bc[c]["HminusD"]["mean"] for c in rr]
    cols = [RED if bc[c]["is_lepr"] else GREY for c in rr]
    ax[0].barh(yv, mean, xerr=[lo, hi], color=cols, edgecolor="white", lw=0.3,
               error_kw={"lw": 0.6, "ecolor": "0.35"})
    ax[0].set_yticks(yv); ax[0].set_yticklabels(rr, fontsize=5.5)
    ax[0].axvline(0, color="k", lw=0.6)
    ax[0].set_xlabel("hyperpolarizing $-$ depolarizing (per-gene SD)")
    ax[0].legend(handles=[Patch(color=RED, label="LepR cluster"),
                          Patch(color=GREY, label="other VMH")],
                 loc="lower right", frameon=False, fontsize=6)
    ax[0].set_title("A   Effector balance across VMH clusters", loc="left")

    pair = cc["lepr_class_vs_others"]
    kk = sorted(pair, key=lambda k: pair[k]["diff_mean_over_seeds"])
    ax[1].barh(np.arange(len(kk)), [pair[k]["diff_mean_over_seeds"] for k in kk],
               xerr=[pair[k]["diff_range"] / 2 for k in kk], capsize=3,
               color=[RED if pair[k]["sign_consistent_across_seeds"] else GREY
                      for k in kk],
               edgecolor="white")
    ax[1].set_yticks(np.arange(len(kk)))
    ax[1].set_yticklabels([k.split("_vs_")[1] for k in kk], fontsize=7)
    ax[1].axvline(0, color="k", lw=.6)
    ax[1].set_xlabel("LepR class $-$ other (H$-$D)")
    ax[1].set_title("B   Class contrasts on the effector axis", loc="left")

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
    ax[0].set_ylabel("Excess over chance (%)")
    ax[0].set_title("A   Within-class structure, by class", loc="left")

    thr = mc["thresholds"]
    for c, col in (("Nfib", BLUE), ("Lepr", RED), ("Esr1", "0.45")):
        v = [mc["by_class"][c][str(t_)].get("excess") for t_ in thr]
        ok = [(t_, x_) for t_, x_ in zip(thr, v) if x_ is not None]
        ax[1].plot([a for a, _ in ok], [b for _, b in ok], "o-", color=col, label=c, ms=4)
    ax[1].axhline(0, color="k", lw=.6)
    ax[1].set_xlabel("Minimum cluster size (nuclei)")
    ax[1].set_ylabel("Excess over chance (%)")
    ax[1].legend(frameon=False, fontsize=7)
    ax[1].set_title("B   Sensitivity to minimum cluster size", loc="left")

    lev = npart["power_matched_control"]["levels"]
    ax[2].plot([l["effect_sd_units"] for l in lev],
               [l["recovered_pct"]["mean"] for l in lev], "o-", color=ORANGE)
    obs = npart["within_lepr_class_nine_clusters"]["pct"]["mean"]
    ax[2].axhline(obs, color=RED, ls="--", lw=.9)
    ax[2].text(0.40, obs, "observed ", color=RED, fontsize=6.5, va="bottom", ha="right")
    ax[2].set_xlabel("Injected cluster effect (SD)")
    ax[2].set_ylabel("Variance recovered (%)")
    ax[2].set_title("C   Detectability given cluster labels", loc="left")

    fig.tight_layout(); save(fig, "Fig3_within_class_survey"); return fig


# ------------------------------------------------------------------ Figure 4
def fig4():
    disc, dep = j("discreteness.json"), j("trpc5_depth.json")
    null = [float(r["lr"]) for r in t("trpc5_depth_null.csv")]
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.3))

    # per-nucleus axis values: discreteness.json stores summaries only, so these come
    # from the npz, which holds ONE control-gene draw. Every other panel here is a
    # 20-draw quantity. The legend says so and so does the axis label.
    _z = np.load(os.path.join(C.DATA, "percell_scores.npz"), allow_pickle=True)
    yv = np.asarray(_z["HminusD"], float)

    ax[0].hist(yv, bins=60, density=True, color="#CFD9E8", edgecolor="white", lw=.3)
    try:
        from scipy.stats import gaussian_kde
        xs = np.linspace(yv.min(), yv.max(), 400)
        ax[0].plot(xs, gaussian_kde(yv)(xs), color=BLUE, lw=1.6, label="observed")
    except Exception:
        pass

    # the smallest split the simulation says would have been recovered, for scale
    # The overlay shows the smallest BALANCED split mixture selection recovers, drawn
    # for scale. It comes from power_grid.json now: 04_discreteness.py no longer
    # writes a grid, and the balanced row is the one this overlay depicts.
    pg = j("power_grid.json")
    sep = pg["thresholds_80pct"]["balanced_2"]["bic"]
    if sep:
        sd = yv.std()
        mu = yv.mean()
        try:
            from scipy.stats import norm
            xs2 = np.linspace(yv.min(), yv.max(), 400)
            comp = 0.5 * norm.pdf(xs2, mu - sep * sd / 2, sd * 0.6) + \
                   0.5 * norm.pdf(xs2, mu + sep * sd / 2, sd * 0.6)
            ax[0].plot(xs2, comp, color=GREY, lw=1.2, ls="--",
                       label=f"{sep:g} SD balanced split (recovered)")
        except Exception:
            pass

    ax[0].axvline(0, color="k", lw=.6)
    d = disc["observed"].get("dip", {})
    if isinstance(d, dict) and "p" in d:
        ax[0].annotate(f"dip D = {d['D']:.4f}, p = {d['p']:.3f}",
                       (0.03, 0.95), xycoords="axes fraction", va="top",
                       fontsize=6.5, color="0.25")
    ax[0].set_xlabel("hyperpolarizing $-$ depolarizing (per-gene SD, one draw)")
    ax[0].set_ylabel("Density")
    ax[0].legend(frameon=False, fontsize=6, loc="upper right")
    ax[0].set_title("A   Per-nucleus effector balance", loc="left")

    # Three configurations, not one curve. The section argues from the CONTRAST
    # between them: the balanced split is recovered, the configuration the
    # subpopulation account predicts is not at any separation, and recovery returns
    # once the displacements are unbalanced enough. A single curve cannot show that.
    seps = pg["separations_sd"]
    for cname, col, lab in (
            ("balanced_2", ORANGE, "50/50 split"),
            ("sohn_3_skew", BLUE, "18/18/64, unbalanced ($-$0.2 / +0.8)"),
            ("sohn_3", RED, "18/18/64, symmetric")):
        y = [pg["grid"][f"{cname}_sd{e}"]["pct_bic"] for e in seps]
        ax[1].plot(seps, y, "s--", color=col, ms=4, label=lab)
    ax[1].axhline(80, color="0.6", ls=":", lw=.7)
    ax[1].set_xlabel("Simulated separation (SD)")
    ax[1].set_ylabel("Simulations recovering two components (%)")
    ax[1].set_ylim(-4, 105)
    ax[1].legend(frameon=False, fontsize=5.5, loc="upper left")
    ax[1].set_title("B   Mixture recovery by configuration", loc="left")

    raw = dep["downsampling"]["pct_detected_raw"]
    ds = dep["downsampling"]["pct_detected_downsampled"]
    ks2 = sorted(raw, key=lambda k: -raw[k]); x = np.arange(len(ks2)); w = .38
    ax[2].bar(x - w/2, [raw[k] for k in ks2], w, label="observed", color=ORANGE)
    ax[2].bar(x + w/2, [ds[k] for k in ks2], w, color=GREY,
              label=f"depth-matched ({dep['downsampling']['target_umi']} UMI)")
    ax[2].set_xticks(x); ax[2].set_xticklabels(ks2, rotation=90, fontsize=6)
    ax[2].set_ylabel(f"Nuclei {it('Trpc5')}-positive (%)")
    lg = dep["logistic_cluster_given_depth"]
    ax[2].legend(frameon=False, fontsize=6.5)
    ax[2].set_title("C   " + r"$\it{Trpc5}$" + " detection by cluster, with depth matching",
                    loc="left")

    print("  note: Fig4A is built from one control-gene draw "
          "(percell_scores.npz); 4B and 4C are 20-draw quantities")
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
    ax[0].set_title("A   Translatome enrichment",
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
    ax[1].set_title("B   Cross-population placement", loc="left")

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


    # ---- record the hash of everything written here, so the record cannot drift
    import hashlib
    produced = []
    for nm in WRITTEN:
        for ext in ("png", "pdf"):
            produced.append(f"{nm}.{ext}")
    produced.append("manuscript_figures.pdf")

    lines = []
    for fn in produced:
        fp = os.path.join(C.FIGURES, fn)
        if os.path.exists(fp):
            h = hashlib.md5(open(fp, "rb").read()).hexdigest()
            lines.append(f"{h}  {fn}")
    with open(os.path.join(C.FIGURES, "checksums.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nchecksums.txt: {len(lines)} files")

    # ---- what else is sitting in figures/ that this script did not write?
    known = set(produced) | {"checksums.txt"}
    extra = sorted(f for f in os.listdir(C.FIGURES)
                   if os.path.isfile(os.path.join(C.FIGURES, f)) and f not in known)
    if extra:
        print(f"\n{len(extra)} file(s) in figures/ not written by this script:")
        for f in extra:
            print(f"    {f}")
        print("  07_figures.py writes the superseded set here too (see README), so this")
        print("  is not necessarily wrong. It is printed because the directory alone")
        print("  does not say which files the manuscript uses, and picking the wrong")
        print("  Fig4 at submission is a mistake nothing else would catch.")


if __name__ == "__main__":
    main()
