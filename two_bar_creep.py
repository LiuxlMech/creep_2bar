#!/usr/bin/env python3
"""Two-bar thermo-creep coupling calculator.

Model:
- Bar A: room temperature, elastic only.
- Bar B: high temperature, elastic + thermal + creep.
- Equal axial strain due to rigid connection, zero resultant force at connection.

Usage examples:
    python3 two_bar_creep.py --hours 500 --k 2
    python3 two_bar_creep.py --hours 500 --k 2 --creep-model creep_models/garofalo.py
"""

from __future__ import annotations

import argparse
import importlib.util
import math
from pathlib import Path
from dataclasses import dataclass
from types import ModuleType
from typing import Callable, Dict, List, Tuple


# (Temperature °C, mean alpha from 20°C to T, unit 1/°C)
ALPHA_DATA: List[Tuple[float, float]] = [
    (204.0, 1.31e-5),
    (316.0, 1.35e-5),
    (427.0, 1.39e-5),
    (538.0, 1.44e-5),
    (649.0, 1.51e-5),
    (760.0, 1.57e-5),
    (871.0, 1.66e-5),
    (982.0, 1.73e-5),
]

# (Temperature °C, E MPa)
E_DATA: List[Tuple[float, float]] = [
    (20.0, 208000.0),
    (93.0, 201000.0),
    (204.0, 199000.0),
    (316.0, 191000.0),
    (427.0, 184000.0),
    (538.0, 176000.0),
    (649.0, 168000.0),
    (760.0, 157000.0),
    (871.0, 146000.0),
    (982.0, 129000.0),
]


def linear_interp(x: float, table: List[Tuple[float, float]]) -> float:
    """Piecewise linear interpolation with endpoint clamping."""
    if x <= table[0][0]:
        return table[0][1]
    if x >= table[-1][0]:
        return table[-1][1]
    for i in range(len(table) - 1):
        x0, y0 = table[i]
        x1, y1 = table[i + 1]
        if x0 <= x <= x1:
            r = (x - x0) / (x1 - x0)
            return y0 + r * (y1 - y0)
    raise RuntimeError("Interpolation failed unexpectedly.")


def alpha_mean_from_20_at(temp_c: float) -> float:
    """Mean thermal expansion coefficient from 20°C to temp_c."""
    return linear_interp(temp_c, ALPHA_DATA)


def E_at(temp_c: float) -> float:
    return linear_interp(temp_c, E_DATA)


def _build_cumulative_thermal_strain_table() -> List[Tuple[float, float]]:
    """Build (T, eps_th_from_20_to_T) from mean-alpha data."""
    table: List[Tuple[float, float]] = [(20.0, 0.0)]
    for t, alpha_mean in ALPHA_DATA:
        eps = alpha_mean * (t - 20.0)
        table.append((t, eps))
    return table


EPS_FROM_20_TABLE = _build_cumulative_thermal_strain_table()


def thermal_strain_from_20(temp_c: float) -> float:
    """Thermal strain from 20°C to temp_c, reconstructed from mean-alpha table."""
    if temp_c <= 20.0:
        return 0.0
    # For 20~204°C, use alpha_mean(20->204) as constant mean slope near the lower bound.
    if temp_c <= 204.0:
        return ALPHA_DATA[0][1] * (temp_c - 20.0)
    return linear_interp(temp_c, EPS_FROM_20_TABLE)


def thermal_strain_between(T0_c: float, T1_c: float) -> float:
    """Thermal strain between two temperatures using mean-alpha interpretation."""
    return thermal_strain_from_20(T1_c) - thermal_strain_from_20(T0_c)


RateFunc = Callable[..., float]


@dataclass
class CreepModel:
    name: str
    params: Dict[str, float]
    rate_func: RateFunc
    source_path: str


@dataclass
class SimConfig:
    T_a_c: float = 20.0
    T_b_c: float = 800.0
    hours: float = 500.0
    dt_h: float = 0.01
    epsilon_creep0: float = 0.0
    k_area_ratio: float = 2.0  # k = A_A / A_B
    creep_model_path: str = "creep_models/norton.py"


@dataclass
class SimResult:
    time_h: List[float]
    eps_creep: List[float]
    sigma_b_mpa: List[float]
    sigma_a_mpa: List[float]
    eps_th_b: float
    E_a_mpa: float
    E_b_mpa: float
    k_area_ratio: float
    creep_model_name: str
    creep_model_path: str
    creep_model_params: Dict[str, float]


def _load_python_module(path: Path) -> ModuleType:
    module_name = f"creep_model_{path.stem}_{abs(hash(path.resolve()))}"
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ValueError(f"Failed to load creep model file: {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_creep_model(path_str: str) -> CreepModel:
    path = Path(path_str).expanduser()
    if not path.is_absolute():
        path = Path.cwd() / path
    if not path.exists():
        raise FileNotFoundError(f"Creep model file not found: {path}")
    mod = _load_python_module(path)
    rate_func = getattr(mod, "creep_rate", None)
    if not callable(rate_func):
        raise ValueError(f"{path} must define callable `creep_rate(...)`")
    raw_params = getattr(mod, "PARAMS", {})
    if raw_params is None:
        raw_params = {}
    if not isinstance(raw_params, dict):
        raise ValueError(f"{path} field `PARAMS` must be a dict")
    name = str(getattr(mod, "MODEL_NAME", path.stem))
    params: Dict[str, float] = {str(k): float(v) for k, v in raw_params.items()}
    return CreepModel(name=name, params=params, rate_func=rate_func, source_path=str(path))


class TwoBarCreepSolver:
    def __init__(self, cfg: SimConfig):
        self.cfg = cfg
        self.E_a = E_at(cfg.T_a_c)
        self.E_b = E_at(cfg.T_b_c)
        self.eps_th_b = thermal_strain_between(cfg.T_a_c, cfg.T_b_c)
        self.creep_model = load_creep_model(cfg.creep_model_path)
        if cfg.k_area_ratio <= 0:
            raise ValueError("k_area_ratio must be > 0")
        # Force equilibrium: sigma_A * A_A + sigma_B * A_B = 0, k = A_A/A_B.
        # Compatibility gives sigma_B = - (eps_th + eps_creep) / (1/(k*E_A) + 1/E_B).
        self.C = 1.0 / (1.0 / (cfg.k_area_ratio * self.E_a) + 1.0 / self.E_b)

    def sigma_b(self, eps_creep_b: float) -> float:
        # Compression in B is negative.
        return -self.C * (eps_creep_b + self.eps_th_b)

    def creep_rate(self, eps_creep_b: float, t_h: float) -> float:
        sigma = self.sigma_b(eps_creep_b)
        try:
            value = self.creep_model.rate_func(
                sigma_mpa=sigma,
                eps_creep=eps_creep_b,
                t_h=t_h,
                temp_c=self.cfg.T_b_c,
                params=self.creep_model.params,
            )
        except Exception as exc:
            raise ValueError(
                f"Error in creep model '{self.creep_model.name}' ({self.creep_model.source_path})"
            ) from exc
        if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            raise ValueError("creep_rate(...) must return a finite float")
        return float(value)

    def solve(self) -> SimResult:
        n = int(round(self.cfg.hours / self.cfg.dt_h))
        t = 0.0
        eps_c = self.cfg.epsilon_creep0

        time_h = [t]
        eps_creep = [eps_c]
        sigma_b = [self.sigma_b(eps_c)]
        sigma_a = [-sigma_b[0] / self.cfg.k_area_ratio]

        # RK4 integration for better stability than forward Euler.
        for _ in range(n):
            dt = self.cfg.dt_h

            k1 = self.creep_rate(eps_c, t)
            k2 = self.creep_rate(eps_c + 0.5 * dt * k1, t + 0.5 * dt)
            k3 = self.creep_rate(eps_c + 0.5 * dt * k2, t + 0.5 * dt)
            k4 = self.creep_rate(eps_c + dt * k3, t + dt)
            eps_c = eps_c + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

            t += dt
            sb = self.sigma_b(eps_c)

            time_h.append(t)
            eps_creep.append(eps_c)
            sigma_b.append(sb)
            sigma_a.append(-sb / self.cfg.k_area_ratio)

        return SimResult(
            time_h=time_h,
            eps_creep=eps_creep,
            sigma_b_mpa=sigma_b,
            sigma_a_mpa=sigma_a,
            eps_th_b=self.eps_th_b,
            E_a_mpa=self.E_a,
            E_b_mpa=self.E_b,
            k_area_ratio=self.cfg.k_area_ratio,
            creep_model_name=self.creep_model.name,
            creep_model_path=self.creep_model.source_path,
            creep_model_params=self.creep_model.params,
        )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Two-bar thermo-creep coupling calculator")
    p.add_argument("--T-a", type=float, default=20.0, help="Bar A temperature (°C)")
    p.add_argument("--T-b", type=float, default=800.0, help="Bar B temperature (°C)")
    p.add_argument("--hours", type=float, default=500.0, help="Simulation duration (h)")
    p.add_argument("--k", type=float, default=2.0, help="Area ratio k=A_A/A_B (e.g. 2)")
    p.add_argument("--dt", type=float, default=0.01, help="Time step (h), suggest <=0.1 for stability")
    p.add_argument(
        "--creep-model",
        type=str,
        default="creep_models/norton.py",
        help="Path to creep model file (.py) defining creep_rate(...) and PARAMS",
    )
    return p


def main() -> None:
    args = build_parser().parse_args()
    cfg = SimConfig(
        T_a_c=args.T_a,
        T_b_c=args.T_b,
        hours=args.hours,
        k_area_ratio=args.k,
        dt_h=args.dt,
        creep_model_path=args.creep_model,
    )
    result = TwoBarCreepSolver(cfg).solve()

    print("=== Two-Bar Creep Result ===")
    print(f"E_A = {result.E_a_mpa:.2f} MPa")
    print(f"E_B = {result.E_b_mpa:.2f} MPa")
    print(f"k = A_A/A_B = {result.k_area_ratio:.6f}")
    print(f"creep_model = {result.creep_model_name}")
    print(f"creep_model_file = {result.creep_model_path}")
    if result.creep_model_params:
        print(f"creep_model_params = {result.creep_model_params}")
    print(f"eps_th_B = {result.eps_th_b:.6e}")
    print(f"sigma_B(0) = {result.sigma_b_mpa[0]:.3f} MPa")
    print(f"sigma_B({args.hours:g} h) = {result.sigma_b_mpa[-1]:.3f} MPa")
    print(f"eps_creep_B({args.hours:g} h) = {result.eps_creep[-1]:.6e}")
    if not math.isfinite(result.sigma_b_mpa[-1]) or not math.isfinite(result.eps_creep[-1]):
        print("[WARN] Numerical overflow detected. Please reduce --dt (e.g., 0.01 or 0.001).")


if __name__ == "__main__":
    main()
