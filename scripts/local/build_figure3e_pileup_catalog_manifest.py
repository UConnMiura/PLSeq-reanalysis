#!/usr/bin/env python3
"""Build the 28-event configuration used by the pileup-catalog Slurm array."""

from pathlib import Path
import csv
from collections import defaultdict

ROOT = Path("/Users/pedromiura/Documents/codex/PMID37679364_reanalysis")
REF = ROOT / "reference/ref"
ATLAS = ROOT / "figures/figure3e_loci/atlas/figure3e_locus_atlas_manifest.tsv"
COUNTS = ROOT / "hpc_results/figure3e/figure3e_reconstructed_counts.tsv"
OUT = ROOT / "metadata/figure3e_pileup_catalog_manifest.tsv"


def bed(path):
    ans = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            ans[f[3]] = (f[0], int(f[1]), int(f[2]))
    return ans


def main():
    beds = {
        "up": bed(REF / "up.figure3e.reconstructed.bed"),
        "down": bed(REF / "down.figure3e.reconstructed.bed"),
        "common": bed(REF / "UTR.figure3e.bed"),
        "distal": bed(REF / "dis.figure3e.bed"),
        "ce": bed(REF / "SE_PSI.figure3e.suppdata9.bed"),
    }

    event_meta = {}
    with ATLAS.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            event_meta.setdefault(r["gene"], (int(r["order"]), r["bed_target"], r["strand"]))

    depth = defaultdict(int)
    with COUNTS.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["isoform"] in {"Long", "Short"}:
                depth[(r["gene"], r["probe_set"])] += int(r["total"])

    fields = [
        "array_index", "order", "gene", "event", "probe_set", "strand", "chrom",
        "view_start", "view_end", "ce_start", "ce_end",
        "up_start", "up_end", "down_start", "down_end",
        "common_start", "common_end", "distal_start", "distal_end",
        "selection_depth",
    ]
    rows = []
    for gene, (order, event, strand) in event_meta.items():
        probe = max(("p1", "p2"), key=lambda p: depth[(gene, p)])
        chroms = {beds[k][event if k in {"up", "down", "ce"} else gene][0] for k in beds}
        if len(chroms) != 1:
            raise ValueError(f"Chromosome mismatch for {gene}: {chroms}")
        chrom = chroms.pop()
        up = beds["up"][event]; down = beds["down"][event]
        common = beds["common"][gene]; distal = beds["distal"][gene]; ce = beds["ce"][event]
        starts = [x[1] for x in (up, down, common, distal, ce)]
        ends = [x[2] for x in (up, down, common, distal, ce)]
        pad = max(75, round((max(ends) - min(starts)) * 0.015))
        rows.append({
            "order": order, "gene": gene, "event": event, "probe_set": probe,
            "strand": strand, "chrom": chrom,
            "view_start": max(0, min(starts) - pad), "view_end": max(ends) + pad,
            "ce_start": ce[1], "ce_end": ce[2],
            "up_start": up[1], "up_end": up[2],
            "down_start": down[1], "down_end": down[2],
            "common_start": common[1], "common_end": common[2],
            "distal_start": distal[1], "distal_end": distal[2],
            "selection_depth": depth[(gene, probe)],
        })
    rows.sort(key=lambda r: r["order"])
    for i, row in enumerate(rows):
        row["array_index"] = i
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"Wrote {len(rows)} events to {OUT}")


if __name__ == "__main__":
    main()
