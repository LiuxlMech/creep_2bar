"""Garofalo hyperbolic-sine creep model example.

Rate law:
    eps_dot = A * sign(sigma) * sinh(alpha * |sigma|)^n
"""

import math

MODEL_NAME = "Garofalo"
PARAMS = {
    "A": 1.0e-9,    # 1/h
    "alpha": 0.001, # 1/MPa
    "n": 3.0,
}


def creep_rate(*, sigma_mpa: float, eps_creep: float, t_h: float, temp_c: float, params: dict) -> float:
    del eps_creep, t_h, temp_c
    A = float(params["A"])
    alpha = float(params["alpha"])
    n = float(params["n"])
    sgn = 1.0 if sigma_mpa >= 0.0 else -1.0
    x = alpha * abs(sigma_mpa)
    # Clamp argument to avoid overflow in example model.
    x = min(x, 50.0)
    return A * sgn * (math.sinh(x) ** n)
