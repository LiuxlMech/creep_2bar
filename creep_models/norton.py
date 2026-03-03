"""Norton creep model.

Interface required by two_bar_creep.py:
- MODEL_NAME: str
- PARAMS: dict[str, float]
- creep_rate(...): returns d(eps_creep)/dt in 1/h
"""

MODEL_NAME = "Norton"
PARAMS = {
    "K": 3.7e-11,  # 1/(h*MPa^n)
    "n": 3.0,
}


def creep_rate(*, sigma_mpa: float, eps_creep: float, t_h: float, temp_c: float, params: dict) -> float:
    del eps_creep, t_h, temp_c
    K = float(params["K"])
    n = float(params["n"])
    # Sign-compatible form: eps_dot = K * sigma * |sigma|^(n-1).
    return K * sigma_mpa * (abs(sigma_mpa) ** (n - 1.0))
