"""Strict two-column readers. Raw files are immutable; roundoff handling is explicit."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
import hashlib
import json
import numpy as np

ZERO_TOLERANCE = 1e-12

@dataclass(frozen=True)
class Curve:
    strain: np.ndarray
    stress: np.ndarray
    source: str = "in_memory"
    sha256: str | None = None
    origin_adjustments: tuple = field(default_factory=tuple)

    def validate(self) -> None:
        x, y = np.asarray(self.strain), np.asarray(self.stress)
        if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) < 3:
            raise ValueError("At least three paired one-dimensional samples are required.")
        if not np.isfinite(x).all() or not np.isfinite(y).all():
            raise ValueError("Non-finite values are not accepted.")
        if np.any(x < 0) or np.any(x >= 1):
            raise ValueError("Expected compression-positive engineering strain in [0, 1).")
        if np.any(np.diff(x) <= 0):
            raise ValueError("Strain must be strictly increasing; loading branches are not sorted or merged.")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def load_curve(path: Path, *, strain_unit: str = "fraction", stress_unit: str = "GPa") -> Curve:
    """Read two comma/whitespace columns; retain negative stress and all loading data.

    Only strain magnitudes < 1e-12 are snapped to zero, recorded in metadata.
    No zero-stress correction, smoothing, sign censoring or extrapolation is done.
    """
    path = Path(path)
    scales = {"GPa": 1., "MPa": 1e-3, "Pa": 1e-9, "bar": 1e-4}
    if strain_unit not in {"fraction", "percent"} or stress_unit not in scales:
        raise ValueError("Specify strain as fraction/percent and stress as GPa/MPa/Pa/bar.")
    rows = []
    header_seen = False
    for line_no, raw in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        tokens = [t.strip() for t in line.split(",")] if "," in line else line.split()
        if len(tokens) != 2:
            raise ValueError(f"{path.name}:{line_no}: expected two columns.")
        if tuple(t.lower() for t in tokens) in {("strain", "stress"), ("strain_fraction", "stress_gpa"), ("strain", "stress_gpa")}:
            if rows or header_seen:
                raise ValueError(f"{path.name}:{line_no}: misplaced/repeated header.")
            header_seen = True
            continue
        try:
            rows.append(tuple(float(t) for t in tokens))
        except ValueError as e:
            raise ValueError(f"{path.name}:{line_no}: invalid numerical record.") from e
    if len(rows) < 3:
        raise ValueError("At least three samples are required.")
    a = np.array(rows, dtype=float)
    x = a[:, 0] / (100. if strain_unit == "percent" else 1.)
    y = a[:, 1] * scales[stress_unit]
    small = (np.abs(x) < ZERO_TOLERANCE) & (x != 0)
    changes = tuple({"sample": int(i), "strain_before": float(x[i]), "strain_after": 0.} for i in np.where(small)[0])
    x[small] = 0.
    curve = Curve(x, y, path.name, sha256_file(path), changes)
    curve.validate()
    return curve
