"""Render the pinned WASP-39 b spectrum and its marginal uncertainties."""

import json
from pathlib import Path

import matplotlib.pyplot as plt


def main():
    root = Path(__file__).resolve().parents[1]
    data = json.loads((root / "data/wasp39b_benchmark.json").read_text())
    colours = {
        "NIRISS_SOSS_Order1": "#005677",
        "NIRISS_SOSS_Order2": "#62428a",
        "NIRCam_F322W2": "#934626",
        "NIRSpec_G395H_NRS1": "#36723a",
        "NIRSpec_G395H_NRS2": "#ad2859",
        "NIRSpec_PRISM": "#575a5f",
    }
    ordered = ["NIRSpec_PRISM", *(name for name in colours if name != "NIRSpec_PRISM")]
    with plt.rc_context({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}):
        fig, (spectrum, errors) = plt.subplots(
            2, 1, sharex=True, height_ratios=(2, 1), figsize=(8.4, 6.0), layout="constrained",
        )
        fig.suptitle("WASP-39 b · JWST transmission spectrum", x=0.01, ha="left", fontweight="bold")
        for name in ordered:
            points = [point for point in data["points"] if point["instrument"] == name]
            wavelength = [point["wavelength_um"] for point in points]
            depth = [point["transit_depth_ppm"] for point in points]
            sigma = [point["uncertainty_ppm"] for point in points]
            label = name.replace("_", " ")
            spectrum.errorbar(
                wavelength, depth, yerr=sigma, fmt="o", ms=2.6,
                color=colours[name], alpha=0.75, lw=0.65, label=f"{label} ({len(points)})",
            )
            errors.scatter(wavelength, sigma, s=9, color=colours[name], alpha=0.75)
        spectrum.set_ylabel("Transit depth (ppm)")
        spectrum.legend(ncol=2, frameon=False, loc="upper left", fontsize=8)
        errors.set(xlabel="Wavelength (µm)", ylabel="Marginal σ (ppm)", yscale="log")
        for axis in (spectrum, errors):
            axis.grid(axis="y", color="#e0e4e8", linewidth=0.65)
            axis.set_axisbelow(True)
        output = root / "assets/wasp39b-transmission-evidence.png"
        output.parent.mkdir(exist_ok=True)
        fig.savefig(output, dpi=115)
        plt.close(fig)


if __name__ == "__main__":
    main()
