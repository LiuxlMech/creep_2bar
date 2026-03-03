"""Time-hardened power-law creep model example.

Rate law:
    eps_dot = A * sign(sigma) * |sigma|^n * (t + t0)^m
"""

MODEL_NAME = "Time-Hardened Power Law"
PARAMS = {
    "A": 1.0e-14,
    "n": 3.0,
    "m": -0.25,
    "t0_h": 1.0,
}


def creep_rate(*, sigma_mpa: float, eps_creep: float, t_h: float, temp_c: float, params: dict) -> float:
    del eps_creep, temp_c
    A = float(params["A"])
    n = float(params["n"])
    m = float(params["m"])
    t0_h = float(params["t0_h"])
    sgn = 1.0 if sigma_mpa >= 0.0 else -1.0
    return A * sgn * (abs(sigma_mpa) ** n) * ((max(t_h, 0.0) + t0_h) ** m)
