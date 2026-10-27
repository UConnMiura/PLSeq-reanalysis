#!/usr/bin/env python3
"""Reanalyse Figure 6C as a genotype-by-3'UTR interaction."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


CE = {
    "Dys": 33, "Khc-73": 15, "pod1": 9, "X11L": 8, "spen": 18,
    "CASK": 17, "Mbs": 11, "mtd": 12, "nuf": 9, "CG8671": 10,
    "tutl": 15, "Camta": 9, "par-1": 13, "Evi5": 2, "shi": 10,
    "CdGAPr": 4, "Ptp10D": 12, "Nc73EF": 15, "gish": 14, "Vav": 12,
    "stai": 6, "CG15312": 6, "Crag": 10, "lap": 15, "Eip63E": 4,
    "Calx": 4, "Lar": 16, "Dscam1": 23,
}

CONTROL = "#5B8FA8"
ELAV = "#D89B35"
TEAL = "#178C8C"
PURPLE = "#7A5195"
GRAY = "#73808C"
SIGNIFICANT = "#E45756"


def normalize_sample(value: str) -> str:
    value = value.strip().lower()
    return "Control" if value == "control" else "elav5"


def read_left(source: Path, sheet: int, value_col: int, sample_col: int):
    path = source / f"41467-2023-41207-moesm15-esm-sheet-{sheet}.csv"
    result = defaultdict(list)
    with path.open(newline="") as handle:
        rows = csv.reader(handle)
        next(rows)
        for row in rows:
            try:
                gene = row[0].strip()
                # Supplementary sheet 12 labels the four mutant Calx rows
                # ``Calx_2``; their values exactly match the Calx conditional
                # long/short differences in sheets 10 and 11.
                if sheet == 12 and gene == "Calx_2":
                    gene = "Calx"
                sample = normalize_sample(row[sample_col])
                value = float(row[value_col])
            except (ValueError, IndexError):
                continue
            if gene:
                result[(gene, sample)].append(value)
    return result


def best_pairs(long_values, short_values, target_abs):
    """Recover replicate pairs by matching the deposited absolute differences."""
    n = len(target_abs)
    if n == 0:
        return None
    best = None
    for li in itertools.combinations(range(len(long_values)), n):
        selected_long = [long_values[i] for i in li]
        for si in itertools.combinations(range(len(short_values)), n):
            selected_short = [short_values[i] for i in si]
            for perm in itertools.permutations(selected_short):
                signed = np.array(selected_long) - np.array(perm)
                error = max(abs(a-b) for a, b in zip(sorted(abs(signed)), sorted(target_abs)))
                candidate = (error, selected_long, list(perm), signed)
                if best is None or candidate[0] < best[0]:
                    best = candidate
    if best is None or best[0] > 2e-6:
        return None
    return best


def welch(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2 or len(y) < 2:
        return np.nan, np.nan, np.nan, np.nan
    estimate = np.mean(y) - np.mean(x)
    test = stats.ttest_ind(y, x, equal_var=False)
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    se = math.sqrt(vx/len(x) + vy/len(y))
    df = (vx/len(x) + vy/len(y))**2 / (
        (vx/len(x))**2/(len(x)-1) + (vy/len(y))**2/(len(y)-1)
    ) if se else np.inf
    crit = stats.t.ppf(0.975, df) if np.isfinite(df) else 1.96
    return estimate, estimate-crit*se, estimate+crit*se, float(test.pvalue)


def bh(pvalues):
    p = np.asarray(pvalues, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[ok])]
    ranked = p[order] * len(order) / np.arange(1, len(order)+1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.minimum(ranked, 1)
    return out


def qlabel(q):
    return "NA" if not np.isfinite(q) else (f"{q:.1e}" if q < .001 else f"{q:.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    long = read_left(args.source, 10, 3, 2)
    short = read_left(args.source, 11, 3, 2)
    absdiff = read_left(args.source, 12, 1, 2)
    dpau = read_left(args.source, 8, 2, 1)
    psiall = read_left(args.source, 9, 3, 2)

    # Restrict analysis to events represented in both conditional-PSI source
    # tables and both genotypes.
    candidate_genes = {g for g, _ in absdiff}
    genes = sorted(
        gene for gene in candidate_genes
        if all(long.get((gene, sample)) and short.get((gene, sample))
               for sample in ("Control", "elav5"))
    )
    excluded = sorted(candidate_genes - set(genes))
    summaries, pair_rows = [], []
    for gene in genes:
        paired = {}
        for sample in ("Control", "elav5"):
            match = best_pairs(long[(gene, sample)], short[(gene, sample)],
                               absdiff[(gene, sample)])
            if match is None:
                continue
            error, lv, sv, signed = match
            paired[sample] = signed
            for i, (lval, sval, dval) in enumerate(zip(lv, sv, signed), 1):
                pair_rows.append({
                    "gene": gene, "CE": CE.get(gene, "NA"), "sample": sample,
                    "reconstructed_replicate": i, "PSI_long": lval,
                    "PSI_short": sval, "signed_coupling_long_minus_short": dval,
                    "absolute_coupling": abs(dval), "matching_error": error,
                })
        cmean = np.mean(long[(gene, "Control")]) - np.mean(short[(gene, "Control")])
        emean = np.mean(long[(gene, "elav5")]) - np.mean(short[(gene, "elav5")])
        if set(paired) == {"Control", "elav5"}:
            interaction, ilo, ihi, ip = welch(paired["Control"], paired["elav5"])
        else:
            # Retain incomplete genes descriptively, without fabricating a paired variance.
            interaction, ilo, ihi, ip = emean - cmean, np.nan, np.nan, np.nan
        abs_change, alo, ahi, apval = welch(absdiff[(gene, "Control")],
                                            absdiff[(gene, "elav5")])
        dp_change, dp_lo, dp_hi, dp_p = welch(
            dpau[(gene, "Control")], dpau[(gene, "elav5")])
        all_change, all_lo, all_hi, all_p = welch(
            psiall[(gene, "Control")], psiall[(gene, "elav5")])
        if np.sign(cmean) != np.sign(emean) and abs(cmean) > .02 and abs(emean) > .02:
            pattern = "reversal"
        elif abs(emean) < abs(cmean):
            pattern = "convergence"
        else:
            pattern = "amplification"
        summaries.append({
            "gene": gene, "CE": CE.get(gene, "NA"),
            "n_control": len(absdiff[(gene, "Control")]),
            "n_elav5": len(absdiff[(gene, "elav5")]),
            "control_mean_short": np.mean(short[(gene, "Control")]),
            "control_mean_long": np.mean(long[(gene, "Control")]),
            "elav5_mean_short": np.mean(short[(gene, "elav5")]),
            "elav5_mean_long": np.mean(long[(gene, "elav5")]),
            "control_signed_coupling": cmean, "elav5_signed_coupling": emean,
            "interaction_elav5_minus_control": interaction,
            "interaction_ci_low": ilo, "interaction_ci_high": ihi,
            "interaction_p": ip,
            "absolute_coupling_change": abs_change,
            "absolute_change_ci_low": alo, "absolute_change_ci_high": ahi,
            "absolute_change_p": apval,
            "dPAU_change": dp_change, "dPAU_ci_low": dp_lo,
            "dPAU_ci_high": dp_hi, "dPAU_p": dp_p,
            "PSIAll_change": all_change, "PSIAll_ci_low": all_lo,
            "PSIAll_ci_high": all_hi, "PSIAll_p": all_p,
            "pattern": pattern,
        })

    for field in ("interaction", "absolute_change", "dPAU", "PSIAll"):
        qs = bh([row[f"{field}_p"] for row in summaries])
        for row, q in zip(summaries, qs):
            row[f"{field}_BH_q"] = q

    summaries.sort(key=lambda r: r["absolute_coupling_change"])

    summary_path = args.output / "figure6c_interaction_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(summaries[0]))
        writer.writeheader(); writer.writerows(summaries)
    with (args.output / "figure6c_reconstructed_pairs.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(pair_rows[0]))
        writer.writeheader(); writer.writerows(pair_rows)

    # New Figure 6C: conditional PSI geometry, interaction forest, and marginal context.
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    n = len(summaries)
    fig = plt.figure(figsize=(17.2, 11.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[2.7, 2.0, 2.75], wspace=.34)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    cgs = gs[0, 2].subgridspec(1, 2, wspace=.55)
    ax3 = fig.add_subplot(cgs[0, 0])
    ax4 = fig.add_subplot(cgs[0, 1])
    y = np.arange(n)

    for i, row in enumerate(summaries):
        # Offset genotype rows; circles are short and triangles are long.
        for offset, prefix, color in ((-.13, "control", CONTROL), (.13, "elav5", ELAV)):
            s = row[f"{prefix}_mean_short"]; l = row[f"{prefix}_mean_long"]
            ax1.plot([s, l], [i+offset, i+offset], color=color, lw=1.6, alpha=.9)
            ax1.scatter(s, i+offset, s=24, marker="o", color=color, zorder=3)
            ax1.scatter(l, i+offset, s=28, marker="^", color=color, zorder=3)
    ax1.set_xlim(-.03, 1.03); ax1.set_ylim(n-.5, -.5)
    ax1.set_yticks(y, [f"{r['gene']}  CE{r['CE']}" for r in summaries])
    ax1.set_xlabel("Cassette-exon PSI")
    ax1.set_title("A  Conditional PSI within each 3′UTR class", loc="left", weight="bold")
    ax1.grid(axis="x", color="#E4E7EB", lw=.7)
    ax1.scatter([], [], color=CONTROL, label="Control", s=25)
    ax1.scatter([], [], color=ELAV, label=r"$elav^5$", s=25)
    ax1.scatter([], [], color="#555", marker="o", label="Short 3′UTR", s=25)
    ax1.scatter([], [], color="#555", marker="^", label="Long 3′UTR", s=28)
    ax1.legend(ncol=2, frameon=False, loc="lower center", bbox_to_anchor=(.5, 1.01))

    for i, row in enumerate(summaries):
        effect = row["interaction_elav5_minus_control"]
        significant = (np.isfinite(row["interaction_BH_q"])
                       and row["interaction_BH_q"] < .05)
        color = SIGNIFICANT if significant else GRAY
        ax2.plot([row["interaction_ci_low"], row["interaction_ci_high"]], [i, i],
                 color=color, lw=1.5)
        ax2.scatter(effect, i, color=color, s=28, zorder=3)
        ax2.text(1.02, i, qlabel(row["interaction_BH_q"]),
                 transform=ax2.get_yaxis_transform(), va="center", fontsize=7)
    ax2.axvline(0, color="#87909A", lw=1)
    ax2.set_xlim(-1.05, 1.05); ax2.set_ylim(n-.5, -.5); ax2.set_yticks([])
    ax2.set_xlabel("Interaction effect\n[(Long−Short)$_{elav5}$ − (Long−Short)$_{control}$]")
    ax2.set_title("B  Genotype × 3′UTR interaction", loc="left", weight="bold")
    ax2.text(1.02, 1.012, "BH q", transform=ax2.transAxes,
             fontsize=7, weight="bold", ha="left", va="bottom")
    ax2.scatter([], [], color=SIGNIFICANT, label="BH q < 0.05", s=25)
    ax2.scatter([], [], color=GRAY, label="BH q ≥ 0.05", s=25)
    ax2.legend(frameon=False, loc="lower left", bbox_to_anchor=(-.02, 1.01), ncol=2)

    # Marginal APA and pooled-PSI effects, now displayed as inferential forest plots.
    marginal_specs = (
        (ax3, "dPAU", "C  ΔdPAU", SIGNIFICANT, "s"),
        (ax4, "PSIAll", "D  ΔPSIAll", SIGNIFICANT, "D"),
    )
    for ax, metric, title, sig_color, marker in marginal_specs:
        for i, row in enumerate(summaries):
            q = row[f"{metric}_BH_q"]
            significant = np.isfinite(q) and q < .05
            color = sig_color if significant else GRAY
            ax.plot([row[f"{metric}_ci_low"], row[f"{metric}_ci_high"]],
                    [i, i], color=color, lw=1.35)
            ax.scatter(row[f"{metric}_change"], i, marker=marker, s=24,
                       color=color, zorder=3)
            ax.text(1.03, i, qlabel(q), transform=ax.get_yaxis_transform(),
                    va="center", fontsize=6.5)
        ax.axvline(0, color="#87909A", lw=1)
        ax.set_xlim(-.65, .65); ax.set_ylim(n-.5, -.5); ax.set_yticks([])
        ax.set_xlabel(r"$elav^5$ − control")
        ax.set_title(title, loc="left", weight="bold")
        ax.text(1.03, 1.012, "BH q", transform=ax.transAxes,
                fontsize=7, weight="bold", ha="left", va="bottom")
        ax.grid(axis="x", color="#E4E7EB", lw=.7)
        ax.scatter([], [], marker=marker, color=sig_color,
                   label="q < 0.05", s=24)
        ax.scatter([], [], marker=marker, color=GRAY,
                   label="q ≥ 0.05", s=24)
        ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(-.04, 1.01),
                  ncol=2, columnspacing=.9, handletextpad=.4)

    # Ordering, interval definitions, and multiple-testing details are stated
    # in the manuscript legend, where they remain legible at final size.
    fig.subplots_adjust(top=.90, bottom=.12)
    fig.savefig(args.output / "figure6c_interaction_reanalysis.pdf", bbox_inches="tight")
    fig.savefig(args.output / "figure6c_interaction_reanalysis.png", dpi=300,
                bbox_inches="tight")
    plt.close(fig)

    counts = defaultdict(int)
    for row in summaries: counts[row["pattern"]] += 1
    with (args.output / "README.txt").open("w") as handle:
        handle.write(
            "Figure 6C interaction reanalysis\n"
            "=================================\n\n"
            "Primary estimand: genotype-by-3'UTR interaction =\n"
            "(PSI_long - PSI_short)_elav5 - (PSI_long - PSI_short)_control.\n"
            "Long/short replicate pairs were reconstructed by exact matching to the deposited\n"
            "per-replicate absolute-difference sheet; the matching error is recorded.\n"
            "The four sheet-12 rows labelled Calx_2 were reconciled to Calx because their\n"
            "values exactly match the Calx long/short conditional-PSI differences.\n"
            "Welch tests compare signed coupling values between genotypes. BH correction is\n"
            "applied across genes for each metric family.\n\n"
            + ("Excluded orphan labels lacking conditional PSI data: "
               + ", ".join(excluded) + "\n" if excluded else "")
            + "Pattern counts: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())) + "\n\n"
            "Interpretive narrative\n"
            "------------------------\n"
            "Cassette-exon inclusion was strongly associated with 3'UTR isoform identity in\n"
            "control heads, but this relationship was frequently weakened in elav5 mutants\n"
            "(panel A). A formal genotype-by-3'UTR interaction analysis identified significant\n"
            "changes after BH correction for Calx, Khc-73, Lar, Dscam1, Eip63E, Dys, CG8671,\n"
            "X11L, par-1, and stai (panel B). For all ten events, the absolute PSI separation\n"
            "between long- and short-3'UTR reads decreased in the mutant. Thus, the dominant\n"
            "effect is a loss of conditional coupling between cassette-exon choice and 3'UTR\n"
            "choice, rather than a uniform shift toward exon inclusion or skipping.\n\n"
            "Changes in marginal APA or pooled exon inclusion explain only part of this\n"
            "phenotype. Significant dPAU changes were detected for Calx, Khc-73, Dscam1, and\n"
            "stai (panel C), whereas pooled PSI changed for Dscam1, par-1, Camta, stai, shi,\n"
            "and tutl (panel D). Several strong interaction effects therefore occurred without\n"
            "a significant change in either marginal summary. In particular, Dys, Eip63E,\n"
            "CG8671, X11L, and Lar would be missed by analyses considering only total APA or\n"
            "pooled PSI. Conversely, marginal changes without a significant interaction show\n"
            "that altered APA or exon inclusion does not necessarily imply altered coupling.\n\n"
            "Together, the results support a model in which ELAV helps maintain coordinated\n"
            "selection of cassette exons and alternative 3'UTR isoforms. This analysis\n"
            "establishes association, not a direct molecular mechanism.\n"
        )
    print(f"Wrote {len(summaries)}-gene interaction analysis to {args.output}")


if __name__ == "__main__":
    main()
