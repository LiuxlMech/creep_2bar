"""Template for user-defined creep model.

Required objects:
- MODEL_NAME: str
- PARAMS: dict[str, float]
- creep_rate(...): return d(eps_creep)/dt in 1/h
"""

MODEL_NAME = "My Custom Model"
PARAMS = {
    "p1": 1.0,
    "p2": 2.0,
}


def creep_rate(*, sigma_mpa: float, eps_creep: float, t_h: float, temp_c: float, params: dict) -> float:
    p1 = float(params["p1"])
    p2 = float(params["p2"])
    # Replace with your own equation.
    # Example:
    # return p1 * sigma_mpa + p2 * eps_creep
    return p1 * sigma_mpa * 0.0 + p2 * eps_creep * 0.0 + t_h * 0.0 + temp_c * 0.0
