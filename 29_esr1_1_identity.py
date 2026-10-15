#!/usr/bin/env python3
"""
29_esr1_1_identity.py - is Esr1_1 a subtype of the conserved Esr1 class, or something else?

The per-class survey attributes 15.65% within-class variance in effector composition to the
conserved Esr1 class, and removing Esr1_1 collapses that to 0.38%. Script 28 then found that
Esr1_1 differs from the other Esr1 clusters in roughly 8% of the transcriptome, including a
Satb2-high signature and an adhesion and axon-guidance component, while expressing *less*
Esr1 than its classmates (39% against 80% detection).

That pattern admits two readings. Either Esr1_1 is a genuinely distinct subtype within the
class, which is what the manuscript currently implies, or it is a separate cell type that the
conserved-class assignment placed there, in which case the 15.65% is a between-type difference
mislabelled as within-class and the manuscript's wording is wrong.

Five tests, none of which require new data:

  1. Deposited annotation. What Conserved_cluster values do Esr1_1 nuclei actually carry, and
     do they agree with the other Esr1 clusters? A cluster split across conserved assignments
     is a different thing from one assigned cleanly.

  2. Transcriptome-wide neighbour. Correlate the pseudobulk profile of Esr1_1 against the
     pseudobulk of every conserved class. If its nearest neighbour is not Esr1, the class
     assignment is not supported by overall expression. This is the decisive test.

  3. Doublet burden. Doublets sit between types by construction. If Esr1_1 carries a higher
     doublet-flagged fraction than its classmates, part of its distinctness is technical.

  4. Quality metrics. Total counts and genes detected per nucleus, to rule out Esr1_1 being a
     lower- or higher-quality subset whose apparent difference is a depth artifact. Depth was
     already matched at a ratio of 0.99 in script 28, so this extends that to gene counts.

  5. Region markers. Satb2 is expressed in cortex and amygdala rather than being a canonical
     VMH gene. Detection of Satb2 and of a small set of markers for structures bordering the
     VMH is reported per cluster, so that contamination from an adjacent region can be seen
     if it is there.

Output: results/esr1_1_identity.json
"""
import os, json, warnings
from collections import defaultdict, Counter
import numpy as np
import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "e1id"
FOCUS = "Esr1_1"
# markers for structures bordering or overlapping VMH dissections, plus Satb2 itself.
# Not an exhaustive panel; enough to show gross contamination if present.
REGION = ["Satb2", "Tbr1", "Slc17a7", "Gad1", "Gad2", "Agrp", "Pomc", "Th",
          "Nr5a1", "Esr1", "Fezf1", "Nfib", "Dlk1", "Foxp2", "Lepr"]


def main():
    feats, name2row = V.load_features(PREFIX)
    row2name = {v: k for k, v in name2row.items()}
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        raw = (m.get("Conserved_cluster") or "").strip()
        m["_raw_cons"] = raw
        m["Conserved_class"] = ("" if (not raw or raw.startswith("Doublets"))
                                else raw.split("_")[0])
    cg, _, _ = V.write_col2group(PREFIX, meta, bc2col, "Conserved_class", TAG)
    cs, rs = V.pass_totals(PREFIX, cg, TAG)
    colsum, rowsum = V.read_two_col(cs), V.read_two_col(rs)
    c2c = {}
    for l in open(cg):
        a, b = l.split()
        c2c[int(a)] = b
    cols = sorted(c for c in colsum if c in c2c)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    clu = {c: (meta[col2bc[c]].get("VMH_cluster") or "NA").strip() for c in cols}
    cls = {c: c2c[c] for c in cols}
    out = {}

    # ---------------------------------------------------------------- 1
    print(">> TEST 1: deposited Conserved_cluster assignments", flush=True)
    # every nucleus of the focus VMH cluster, including any excluded as doublets
    all_focus = [b for b, m in meta.items()
                 if (m.get("VMH_cluster") or "").strip() == FOCUS]
    cons = Counter((meta[b]["_raw_cons"] or "unassigned") for b in all_focus)
    print(f"   {FOCUS}: {len(all_focus)} nuclei in the deposited VMH_cluster column")
    for k, v in cons.most_common(8):
        print(f"     {k:<28} {v:5d}  ({100*v/len(all_focus):.1f}%)")
    others = [c for c in cols if clu[c].startswith("Esr1_") and clu[c] != FOCUS]
    cons_o = Counter(meta[col2bc[c]]["_raw_cons"] for c in others)
    print(f"   other Esr1_* clusters: {len(others)} nuclei")
    for k, v in cons_o.most_common(5):
        print(f"     {k:<28} {v:5d}  ({100*v/len(others):.1f}%)")
    out["test1_annotation"] = {
        "focus_n_total": len(all_focus),
        "focus_conserved_assignments": dict(cons.most_common(10)),
        "others_conserved_assignments": dict(cons_o.most_common(10))}

    # ---------------------------------------------------------------- 3
    print("\n>> TEST 3: doublet burden", flush=True)
    foc_d = sum(1 for b in all_focus if meta[b]["_raw_cons"].startswith("Doublets"))
    oth_all = [b for b, m in meta.items()
               if (m.get("VMH_cluster") or "").strip().startswith("Esr1_")
               and (m.get("VMH_cluster") or "").strip() != FOCUS]
    oth_d = sum(1 for b in oth_all if meta[b]["_raw_cons"].startswith("Doublets"))
    print(f"   {FOCUS}: {foc_d}/{len(all_focus)} doublet-flagged "
          f"({100*foc_d/max(len(all_focus),1):.1f}%)")
    print(f"   other Esr1_*: {oth_d}/{len(oth_all)} ({100*oth_d/max(len(oth_all),1):.1f}%)")
    out["test3_doublets"] = {"focus_pct": round(100*foc_d/max(len(all_focus),1), 2),
                             "others_pct": round(100*oth_d/max(len(oth_all),1), 2)}

    # ---------------------------------------------------------------- 2, 4, 5
    print("\n>> TEST 2: transcriptome-wide nearest conserved class", flush=True)
    keep = {r for r, t in rowsum.items() if t > 0}
    ent = V.pass_rows(PREFIX, keep, TAG, "all")
    # pseudobulk: summed counts per group, CP10k, log1p
    groups = defaultdict(list)
    for c in cols:
        groups["FOCUS" if clu[c] == FOCUS else cls[c]].append(c)
    gcols = {g: set(v) for g, v in groups.items()}
    tot = {g: sum(colsum[c] for c in v) for g, v in groups.items()}
    pb = {g: defaultdict(float) for g in groups}
    qc_counts, qc_genes = defaultdict(list), defaultdict(int)
    genes_per_cell = Counter()
    for l in open(ent):
        r, c, v = l.split()
        c, r, v = int(c), int(r), float(v)
        if c not in clu:
            continue
        g = "FOCUS" if clu[c] == FOCUS else cls[c]
        pb[g][r] += v
        genes_per_cell[c] += 1
    common = sorted(keep)
    mats = {}
    for g in pb:
        a = np.array([pb[g].get(r, 0.0) for r in common], float)
        mats[g] = np.log1p(a / max(tot[g], 1) * 1e4)
    focus_vec = mats["FOCUS"]
    sims = {}
    for g, v in mats.items():
        if g == "FOCUS":
            continue
        sims[g] = float(np.corrcoef(focus_vec, v)[0, 1])
    order = sorted(sims, key=lambda k: -sims[k])
    for g in order:
        mark = "  <- assigned class" if g == "Esr1" else ""
        print(f"   {g:<8} r = {sims[g]:.4f}{mark}")
    nearest = order[0]
    out["test2_nearest_class"] = {"correlations": {k: round(v, 4) for k, v in sims.items()},
                                  "nearest": nearest, "assigned": "Esr1",
                                  "assignment_supported": nearest == "Esr1"}
    print(f"   nearest class: {nearest}; assignment "
          f"{'supported' if nearest == 'Esr1' else 'NOT supported'} by overall expression")

    print("\n>> TEST 4: quality metrics", flush=True)
    for g in ("FOCUS", "Esr1"):
        cl = [c for c in groups[g]]
        d = np.array([colsum[c] for c in cl], float)
        ng = np.array([genes_per_cell[c] for c in cl], float)
        label = FOCUS if g == "FOCUS" else "other Esr1 class members"
        print(f"   {label:<26} n={len(cl):4d}  median counts {np.median(d):8.0f}  "
              f"median genes {np.median(ng):6.0f}")
        out.setdefault("test4_quality", {})[g] = {
            "n": len(cl), "median_counts": float(np.median(d)),
            "median_genes": float(np.median(ng))}

    print("\n>> TEST 5: region and identity markers, % of nuclei detected", flush=True)
    want = {name2row[g]: g for g in REGION if g in name2row}
    det = defaultdict(lambda: defaultdict(int))
    for l in open(ent):
        r, c, v = l.split()
        r, c = int(r), int(c)
        if r in want and c in clu and float(v) > 0:
            det[want[r]]["FOCUS" if clu[c] == FOCUS else cls[c]] += 1
    hdr = ["FOCUS", "Esr1", "Lepr", "Nfib", "Fezf1", "Foxp2", "Dlk1"]
    hdr = [h for h in hdr if h in groups]
    print("   " + "gene".ljust(10) + "".join(h[:7].rjust(9) for h in hdr))
    tbl = {}
    for g in REGION:
        if g not in det:
            continue
        row = {h: round(100 * det[g][h] / max(len(groups[h]), 1), 1) for h in hdr}
        tbl[g] = row
        print("   " + g.ljust(10) + "".join(f"{row[h]:8.1f}%" for h in hdr))
    out["test5_markers"] = tbl

    json.dump(out, open(os.path.join(cfg.RESULTS, "esr1_1_identity.json"), "w"), indent=2)
    print("\nWritten to results/esr1_1_identity.json")
    print("\nREADING: if TEST 2 puts the nearest class elsewhere, or TEST 1 shows the focus "
          "cluster split across conserved assignments, or TEST 5 shows non-VMH markers "
          "concentrated in the focus cluster, then the 15.65% is a between-type difference "
          "and the manuscript must say so.")


if __name__ == "__main__":
    main()
