"""Reproduce conditional JWST template and covariance-sensitivity diagnostics."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from spectral_robustness import covariance_sweep, leave_one_product_out, load_benchmark

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, default=ROOT / "data/wasp39b_benchmark.json")
    parser.add_argument("--output", type=Path, default=ROOT / "results/spectral_robustness.json")
    parser.add_argument("--plot", type=Path, default=None)
    args = parser.parse_args()
    points = load_benchmark(args.dataset)
    payload = {
        "dataset": str(args.dataset), "data_type": "published reduced spectrum",
        "caveat": "Descriptive Gaussian features; assumed AR(1) within-product covariance, not a physical retrieval.",
        "covariance_sweep": covariance_sweep(points),
        "product_deletion": leave_one_product_out(points),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if args.plot is not None:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4.5))
        rhos = [x["rho_assumed"] for x in payload["covariance_sweep"]]
        for feature in ("h2o", "co2"):
            vals = [x["conditional_1sigma_ppm"][feature] for x in payload["covariance_sweep"]]
            ax.plot(rhos, vals, marker="o", label=feature.upper())
        ax.set(xlabel="Assumed within-product AR(1) correlation", ylabel="Conditional 1σ (ppm)",
               title="WASP-39 b: model-dependent uncertainty sensitivity")
        ax.legend(); ax.grid(alpha=.25); fig.tight_layout()
        args.plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.plot, dpi=160)
        plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
