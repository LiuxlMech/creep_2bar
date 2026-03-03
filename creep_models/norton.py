"""Temperature-dependent Norton model for Inconel 625.

Rate law:
    eps_dot = K(T) * sigma * |sigma|^(n(T)-1)

Data source (high-temperature alloy Inconel 625):
- Ozgowicz et al., Materials 2023, 16, 7724.
  https://www.mdpi.com/1996-1944/16/24/7724
- Table 3 provides (T, sigma, steady creep rate) data at 650/900/980 C.
- For 650 C, stress exponent n=10.8 is cited from Son et al. in the paper text.
- K(T) and n(T) are inferred from those points in this file.

Notes:
- Creep rates in source table are in %/h; converted here to 1/h.
- Source stress range is about 20-400 MPa. For robustness, stress is clamped to
  this range before evaluating the law.
"""

from __future__ import annotations

import math
from typing import List, Tuple

MODEL_NAME = "Norton (Inconel 625, T-dependent)"

# (Temperature C, K in 1/(h*MPa^n), n)
# 650 C point: n = 10.8 (paper text), K from sample-1 (400 MPa, 9.97e-2 %/h).
# 900/980 C points: n and K inferred from two stress-rate points in Table 3.
_NORTON_TABLE: List[Tuple[float, float, float]] = [
    (650.0, 7.878567350597269e-32, 10.8),
    (900.0, 1.3505374556355225e-13, 4.702779785388754),
    (980.0, 2.678543297487419e-9, 2.4568099394644345),
]

PARAMS = {
    "sigma_min_mpa": 20.0,
    "sigma_max_mpa": 400.0,
}


def _interp_linear_clamped(x: float, points: List[Tuple[float, float]]) -> float:
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    for i in range(len(points) - 1):
        x0, y0 = points[i]
        x1, y1 = points[i + 1]
        if x0 <= x <= x1:
            r = (x - x0) / (x1 - x0)
            return y0 + r * (y1 - y0)
    raise RuntimeError("Interpolation failed unexpectedly.")


def _interp_log_clamped(x: float, points: List[Tuple[float, float]]) -> float:
    log_points = [(tx, math.log(max(ky, 1e-300))) for tx, ky in points]
    return math.exp(_interp_linear_clamped(x, log_points))


def _params_at_temp(temp_c: float) -> tuple[float, float]:
    k_points = [(t, k) for t, k, _ in _NORTON_TABLE]
    n_points = [(t, n) for t, _, n in _NORTON_TABLE]
    K = _interp_log_clamped(temp_c, k_points)
    n = _interp_linear_clamped(temp_c, n_points)
    return K, n


def creep_rate(*, sigma_mpa: float, eps_creep: float, t_h: float, temp_c: float, params: dict) -> float:
    del eps_creep, t_h
    K, n = _params_at_temp(temp_c)
    smin = float(params.get("sigma_min_mpa", 20.0))
    smax = float(params.get("sigma_max_mpa", 400.0))
    s_abs = min(max(abs(sigma_mpa), smin), smax)
    s_eff = s_abs if sigma_mpa >= 0.0 else -s_abs
    # Sign-compatible form: eps_dot = K * sigma * |sigma|^(n-1).
    return K * s_eff * (abs(s_eff) ** (n - 1.0))
