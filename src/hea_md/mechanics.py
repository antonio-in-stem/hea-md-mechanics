"""Finite-window stiffness, operational proof stress, and explicitly bounded maxima."""
from __future__ import annotations
from dataclasses import dataclass, asdict
import math
import numpy as np
from .io import Curve

@dataclass(frozen=True)
class AnalysisConfig:
    fit_min: float = 0.
    fit_max: float = .02
    analysis_max: float = .1
    offset: float = .002
    smooth_points: int = 21
    confirm_points: int = 3

    def validate(self) -> None:
        if not all(math.isfinite(v) for v in [self.fit_min, self.fit_max, self.analysis_max, self.offset]):
            raise ValueError("Configuration values must be finite.")
        if not 0 <= self.fit_min < self.fit_max < self.analysis_max < 1:
            raise ValueError("Require 0 <= fit_min < fit_max < analysis_max < 1.")
        if not 0 < self.offset < self.analysis_max:
            raise ValueError("Offset must be positive and smaller than analysis_max.")
        for name, n in [("smooth_points", self.smooth_points), ("confirm_points", self.confirm_points)]:
            if isinstance(n, bool) or not isinstance(n, int) or n < 1:
                raise ValueError(f"{name} must be a positive integer.")
        if self.smooth_points % 2 != 1:
            raise ValueError("smooth_points must be odd.")


def elastic_fit(curve: Curve, lo: float, hi: float) -> dict:
    curve.validate()
    if not all(math.isfinite(v) for v in [lo, hi]) or not 0 <= lo < hi < 1:
        raise ValueError("Invalid fit interval.")
    if lo < curve.strain[0] - 1e-12 or hi > curve.strain[-1] + 1e-12:
        raise ValueError("Input does not cover the requested fit interval.")
    mask = (curve.strain >= lo - 1e-12) & (curve.strain <= hi + 1e-12)
    x, y = curve.strain[mask], curve.stress[mask]
    if len(x) < 3:
        raise ValueError("Fit requires at least three samples.")
    xc, yc = x - x.mean(), y - y.mean()
    denom = float(xc @ xc)
    if denom <= 0:
        raise ValueError("Fit interval has no strain variation.")
    slope = float((xc @ yc) / denom)
    intercept = float(y.mean() - slope * x.mean())
    if not math.isfinite(slope) or slope <= 0:
        raise ValueError("Fit slope must be positive and finite.")
    r = y - (slope * x + intercept)
    rss, tss = float(r @ r), float(yc @ yc)
    return {"slope_GPa": slope, "intercept_GPa": intercept, "fit_min": lo, "fit_max": hi,
            "n_fit": len(x), "r_squared": None if tss == 0 else 1 - rss / tss,
            "rmse_GPa": float(np.sqrt(rss / len(x))),
            "maximum_absolute_residual_GPa": float(np.max(np.abs(r)))}


def centered_mean(curve: Curve, points: int) -> tuple[np.ndarray, np.ndarray]:
    """Equal-weight centered moving mean, valid samples only; no padding at ends."""
    curve.validate()
    if isinstance(points, bool) or not isinstance(points, int) or points < 1 or points % 2 != 1:
        raise ValueError("Use a positive odd number of averaging samples.")
    if points > len(curve.strain) - 2:
        raise ValueError("Averaging window leaves fewer than three samples.")
    if points == 1:
        return curve.strain.copy(), curve.stress.copy()
    dx = np.diff(curve.strain)
    if not np.allclose(dx, np.median(dx), rtol=1e-7, atol=1e-12):
        raise ValueError("Sample-based averaging requires uniformly spaced strain.")
    half = points // 2
    return curve.strain[half:-half].copy(), np.convolve(curve.stress, np.ones(points)/points, mode="valid")


def offset_crossing(x: np.ndarray, y: np.ndarray, slope: float, intercept: float,
                    *, offset: float, search_min: float, search_max: float, confirm_points: int) -> dict:
    """First persistent downward intersection with E*(strain-offset)+intercept.

    Confirmation suppresses single-sample excursions but does not establish
    irreversible deformation. Every persistent candidate is returned.
    """
    d = y - (slope * (x - offset) + intercept)
    candidates = []
    for i in range(len(x)-1):
        if x[i+1] < search_min - 1e-12 or x[i+1] > search_max + 1e-12:
            continue
        if not (d[i] > 0 and d[i+1] <= 0):
            continue
        end = i+1+confirm_points
        if end > len(x) or x[end-1] > search_max+1e-12 or not np.all(d[i+1:end] <= 0):
            continue
        weight = float(d[i] / (d[i] - d[i+1]))
        cross_x = float(x[i] + weight * (x[i+1]-x[i]))
        if cross_x < search_min - 1e-12:
            continue
        cross_y = float(y[i] + weight * (y[i+1]-y[i]))
        candidates.append({"strain": cross_x, "stress_GPa": cross_y,
                           "bracket": [float(x[i]), float(x[i+1])]})
    return {"status": "found" if candidates else "not_found", "point": candidates[0] if candidates else None,
            "candidate_count": len(candidates), "candidates": candidates}


def peak(x: np.ndarray, y: np.ndarray) -> dict:
    if len(x) == 0:
        raise ValueError("Cannot find peak in an empty window.")
    i = int(np.argmax(y))
    return {"stress_GPa": float(y[i]), "strain": float(x[i]),
            "at_boundary": i in {0, len(x)-1}, "samples": len(x)}


def analyze(curve: Curve, config: AnalysisConfig = AnalysisConfig()) -> dict:
    config.validate()
    curve.validate()
    if curve.strain[0] > 1e-12 or curve.strain[-1] < config.analysis_max - 1e-12:
        raise ValueError("Curve must cover zero strain and the full common comparison window.")
    fit = elastic_fit(curve, config.fit_min, config.fit_max)
    mask = curve.strain <= config.analysis_max + 1e-12
    # Limit all proof-stress processing to the comparison window before averaging.
    window = Curve(curve.strain[mask], curve.stress[mask])
    x, y = centered_mean(window, config.smooth_points)
    if len(x) < config.confirm_points + 2 or x[-1] <= config.fit_max + config.offset:
        raise ValueError("Insufficient search interval after averaging and fitting.")
    e, b = fit["slope_GPa"], fit["intercept_GPa"]
    inside_fit = (x >= config.fit_min) & (x <= config.fit_max)
    contaminated = bool(np.any(y[inside_fit] < e * (x[inside_fit] - config.offset) + b))
    proof = offset_crossing(x,y,e,b,offset=config.offset, search_min=config.fit_max,
                            search_max=config.analysis_max,confirm_points=config.confirm_points)
    if contaminated:
        proof = {"status":"fit_interval_exceeds_offset", "point":None,
                 "candidate_count":proof["candidate_count"],"candidates":proof["candidates"]}
    warnings=[]
    if contaminated: warnings.append("Offset-sized downward excursion inside fit window; no proof stress assigned.")
    if proof["candidate_count"] > 1: warnings.append("Multiple persistent crossings; inspect signal and sensitivity.")
    if proof["point"] is None: warnings.append("No operational proof stress identified under this configuration.")
    raw_peak=peak(window.strain,window.stress)
    if raw_peak['at_boundary']:warnings.append("Common-window maximum is at a boundary, not a resolved interior peak.")
    dx=float(np.median(np.diff(curve.strain)))
    return {"schema_version":1,"input":{"file":curve.source,"sha256":curve.sha256,"samples":len(curve.strain),
                    "negative_stress_samples":int(np.sum(curve.stress<0)),"strain_min":float(curve.strain[0]),"strain_max":float(curve.strain[-1]),
                    "strain_spacing":dx,"origin_roundoff_adjustments":list(curve.origin_adjustments)},
            "config":asdict(config),"elastic_fit":fit,"proof_stress":proof,
            "averaging_span_strain":(config.smooth_points-1)*dx,
            "peak_common_window":raw_peak,"peak_full_record":peak(curve.strain,curve.stress),
            "peak_averaged_common_window":peak(x,y),"warnings":warnings}


def engineering_stress_from_pressure(pzz_bar: float, area_now_A2: float, area_initial_A2: float) -> float:
    if not all(math.isfinite(v) for v in [pzz_bar,area_now_A2,area_initial_A2]) or area_now_A2<=0 or area_initial_A2<=0:
        raise ValueError("Finite pressure and positive areas are required.")
    return 1e-4*pzz_bar*area_now_A2/area_initial_A2
