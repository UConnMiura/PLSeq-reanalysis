#!/usr/bin/env python3
"""Create deterministic PNG and PDF artifacts with the plotting environment."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    x = np.linspace(0, 2 * np.pi, 200)
    fig, ax = plt.subplots(figsize=(4, 2.4))
    ax.plot(x, np.sin(x), color="#5B4B9A", linewidth=2, label="long 3′UTR")
    ax.plot(x, np.cos(x), color="#D18F2E", linewidth=2, label="short 3′UTR")
    ax.set(xlabel="position", ylabel="signal", title="PL-Seq plotting-container smoke test")
    ax.legend(frameon=False)
    fig.tight_layout()

    png = args.output / "plseq_plotting_smoke.png"
    pdf = args.output / "plseq_plotting_smoke.pdf"
    fig.savefig(png, dpi=150)
    with PdfPages(pdf) as document:
        document.savefig(fig)
    plt.close(fig)

    with Image.open(png) as image:
        image.verify()
    if png.stat().st_size < 10_000 or pdf.stat().st_size < 5_000:
        raise SystemExit("Smoke-test artifacts are unexpectedly small")
    print("Wrote {}".format(png))
    print("Wrote {}".format(pdf))


if __name__ == "__main__":
    main()
