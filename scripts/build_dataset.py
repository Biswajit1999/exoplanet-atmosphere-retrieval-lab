"""Build the browser dataset from the published Carter et al. spectrum tables."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "source_products"
OUTPUT = ROOT / "data" / "wasp39b_benchmark.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build() -> dict:
    points: list[dict] = []
    products: list[dict] = []
    excluded_saturated = 0
    for path in sorted(SOURCE.glob("*.csv"), key=lambda item: item.name.casefold()):
        instrument = path.stem
        count = 0
        with path.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                if int(float(row.get("satflag", "0") or 0)) != 0:
                    excluded_saturated += 1
                    continue
                wavelength = float(row["wave"])
                ratio = float(row["rp/rs"])
                ratio_error = 0.5 * (
                    float(row["rp/rs_err_low"]) + float(row["rp/rs_err_hih"])
                )
                points.append(
                    {
                        "wavelength_um": round(wavelength, 8),
                        "wavelength_low_um": round(float(row["wave_low"]), 8),
                        "wavelength_high_um": round(float(row["wave_hig"]), 8),
                        "transit_depth_ppm": round(ratio**2 * 1_000_000, 3),
                        "uncertainty_ppm": round(2 * ratio * ratio_error * 1_000_000, 3),
                        "instrument": instrument,
                    }
                )
                count += 1
        products.append(
            {
                "file": path.name,
                "sha256": sha256(path),
                "included_bins": count,
            }
        )
    points.sort(key=lambda point: (point["wavelength_um"], point["instrument"]))
    return {
        "dataset": "WASP-39 b uniform JWST transmission spectrum",
        "doi": "10.5281/zenodo.10161743",
        "publication_doi": "10.1038/s41550-024-02292-x",
        "archive": {
            "name": "ERS_DataSynthesis_Zenodo.zip",
            "md5": "f617fa28908aa35d334062f1c06bc24d",
            "selection": "Fit_LimbDarkening: R100 for NIRISS, NIRCam and G395H; native bins_scale1 for PRISM",
        },
        "transform": "depth_ppm=(Rp/Rs)^2*1e6; sigma_depth_ppm=2*(Rp/Rs)*mean(asymmetric sigma_RpRs)*1e6",
        "covariance": "No inter-bin covariance matrix is supplied with these selected CSV products; the laboratory uses a diagonal Gaussian likelihood plus an explicit additive error-floor sensitivity parameter.",
        "excluded_saturated_prism_bins": excluded_saturated,
        "products": products,
        "points": points,
    }


if __name__ == "__main__":
    payload = build()
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(f"Wrote {len(payload['points'])} bins to {OUTPUT.relative_to(ROOT)}")
